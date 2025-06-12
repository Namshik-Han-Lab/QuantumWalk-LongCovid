# Potential requirements.txt:
# pandas
# networkx
# pathlib

import itertools
import networkx as nx
import sys
import pickle
import multiprocessing as mp
import pandas as pd
import json
import logging
import math # Added for ceiling division
from pathlib import Path
from typing import List, Tuple, Dict, Any, Iterable, Optional, Set

# --- Configuration ---
# Setup basic logging to provide informative output during execution
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# --- Utility Functions ---
# (save_obj, pairwise, find_shortest_path_edges remain the same as in v2)

def save_obj(obj: Any, file_path: Path) -> None:
    """Saves a Python object to a pickle file."""
    file_addr_pkl = file_path.with_suffix('.pkl')
    try:
        with open(file_addr_pkl, 'wb') as f:
            pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)
        logging.info(f"Successfully saved object to {file_addr_pkl}")
    except (IOError, pickle.PicklingError) as e:
        logging.error(f"Error saving object to {file_addr_pkl}: {e}")
        raise

def pairwise(iterable: Iterable[Any]) -> zip:
    """Generates overlapping pairs from an iterable."""
    a, b = itertools.tee(iterable)
    next(b, None)
    return zip(a, b)

def find_shortest_path_edges(
    pair: Tuple[str, str],
    graph: nx.Graph
) -> List[Tuple[str, str]]:
    """
    Finds all unique edges that constitute *any* shortest path between a source and target node.
    (Worker function for multiprocessing)
    """
    source, target = pair
    all_edges_in_shortest_paths: Set[Tuple[str, str]] = set()

    if source not in graph:
        # logging.warning(f"Source node '{source}' not found in graph within worker for pair {pair}.") # Logged upstream now
        return []
    if target not in graph:
        # logging.warning(f"Target node '{target}' not found in graph within worker for pair {pair}.") # Logged upstream now
        return []
    if source == target:
        return []

    try:
        shortest_paths_iter = nx.all_shortest_paths(graph, source=source, target=target)
        path_found = False
        for path in shortest_paths_iter:
            path_found = True
            for u, v in pairwise(path):
                all_edges_in_shortest_paths.add(tuple(sorted((u, v))))

        if not path_found:
             # logging.debug(f"nx.all_shortest_paths yielded no paths for {pair}, though nodes exist and differ.")
             return []

    except nx.NetworkXNoPath:
        # logging.debug(f"No path found between {source} and {target}.")
        return []
    except nx.NodeNotFound:
        # logging.warning(f"NodeNotFound during shortest path calculation for {pair}. This should have been pre-filtered.")
        return []
    except Exception as e:
        logging.error(f"Unexpected error finding shortest paths for pair {pair}: {e}")
        return []

    return list(all_edges_in_shortest_paths)


# --- Core Logic Functions ---
# (load_config, load_data, build_string_graph, get_protein_lists remain the same as in v2)

def load_config(config_path: str = "config.json") -> Dict[str, Any]:
    """Loads configuration parameters from a JSON file."""
    try:
        config_file = Path(config_path)
        with open(config_file) as f:
            config = json.load(f)
        logging.info(f"Loaded configuration from {config_file}")

        # Set Defaults and Validate
        if "combined_score" not in config:
            config["combined_score"] = 700
            logging.info(f"Using default combined_score: {config['combined_score']}")
        if "num_processes" not in config:
            default_processes = max(1, mp.cpu_count() // 2)
            config["num_processes"] = default_processes
            logging.info(f"Using default num_processes: {config['num_processes']}")
        if "result_dir" not in config:
             config["result_dir"] = "./results"
             logging.info(f"Using default result_dir: {config['result_dir']}")
        # Add chunk_size configuration
        if "chunk_size" not in config:
             config["chunk_size"] = 10000 # Default chunk size for pairs
             logging.info(f"Using default chunk_size: {config['chunk_size']:,}")
        # Add maxtasksperchild configuration
        if "maxtasksperchild" not in config:
             config["maxtasksperchild"] = None # Default: processes live for the life of the pool
             logging.info(f"Using default maxtasksperchild: {config['maxtasksperchild']}")
        else:
             # Ensure it's None or a positive integer
             if config["maxtasksperchild"] is not None and not (isinstance(config["maxtasksperchild"], int) and config["maxtasksperchild"] > 0):
                  logging.warning(f"Invalid maxtasksperchild value ({config['maxtasksperchild']}), setting to None.")
                  config["maxtasksperchild"] = None


        required_keys = ["dip_file", "dep_file", "string_info_file", "string_links_file"]
        missing_keys = [key for key in required_keys if key not in config]
        if missing_keys:
            raise KeyError(f"Missing required key(s) in {config_path}: {', '.join(missing_keys)}")

        config["result_dir"] = Path(config["result_dir"])

        return config
    except FileNotFoundError:
        logging.error(f"Configuration file not found: {config_path}")
        raise
    except json.JSONDecodeError:
        logging.error(f"Error decoding JSON from {config_path}")
        raise
    except KeyError as e:
         logging.error(str(e))
         raise

def load_data(dip_addr: str, dep_addr: str, string_info_addr: str) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Loads the input data files (DIP, DEP, STRING info) into pandas DataFrames."""
    try:
        logging.info(f"Loading DIP data from: {dip_addr}")
        dip_df = pd.read_csv(dip_addr, sep='\t')
        logging.info(f"Loading DEP data from: {dep_addr}")
        dep_df = pd.read_csv(dep_addr, sep='\t')
        logging.info(f"Loading STRING info data from: {string_info_addr}")
        string_info_df = pd.read_csv(string_info_addr, sep='\t')
        logging.info("Successfully loaded DIP, DEP, and STRING info data.")
        return dip_df, dep_df, string_info_df
    except FileNotFoundError as e:
        logging.error(f"Data file not found: {e.filename}")
        raise
    except Exception as e:
        logging.error(f"Error loading data files: {e}")
        raise

def build_string_graph(string_links_addr: str, combined_score_threshold: int) -> nx.Graph:
    """Builds the STRING protein-protein interaction network graph."""
    try:
        logging.info(f"Loading STRING links from: {string_links_addr}")
        string_links_df = pd.read_csv(string_links_addr, sep=' ')
        logging.info(f"Read {len(string_links_df):,} interactions from STRING links file.")

        string_links_filtered = string_links_df[string_links_df['combined_score'] >= combined_score_threshold].copy()
        logging.info(f"Filtered down to {len(string_links_filtered):,} interactions with score >= {combined_score_threshold}")

        if string_links_filtered.empty:
             logging.warning(f"No STRING interactions found with combined_score >= {combined_score_threshold}. Resulting graph will be empty.")
             return nx.Graph()

        string_links_filtered['protein1'] = string_links_filtered['protein1'].astype(str)
        string_links_filtered['protein2'] = string_links_filtered['protein2'].astype(str)

        string_G = nx.from_pandas_edgelist(
            string_links_filtered,
            source='protein1',
            target='protein2',
            edge_attr=['combined_score'],
            create_using=nx.Graph
        )
        logging.info(f"Built STRING graph: {string_G.number_of_nodes():,} nodes, {string_G.number_of_edges():,} edges.")
        return string_G
    except FileNotFoundError:
        logging.error(f"STRING links file not found: {string_links_addr}")
        raise
    except KeyError as e:
        logging.error(f"Missing expected column in {string_links_addr}: {e}")
        raise
    except Exception as e:
        logging.error(f"Error building STRING graph: {e}")
        raise

def get_protein_lists(
    dip_df: pd.DataFrame,
    dep_df: pd.DataFrame,
    string_info_df: pd.DataFrame
) -> Tuple[Dict[str, List[str]], List[str]]:
    """Processes DIP and DEP dataframes to extract relevant STRING protein IDs."""
    logging.info("Mapping DIP/DEP entries to STRING protein IDs...")

    # Process DIP Data
    string_map_df = string_info_df[['preferred_name', '#string_protein_id']].copy()
    string_map_df.rename(columns={'#string_protein_id': 'string_id'}, inplace=True)
    dip_id_df = dip_df.merge(
        string_map_df,
        left_on='StringGene',
        right_on='preferred_name',
        how='inner'
    )
    lost_dips = len(dip_df) - len(dip_id_df)
    if lost_dips > 0:
        logging.warning(f"{lost_dips} DIP entries could not be mapped to STRING IDs via 'StringGene'/'preferred_name' and were excluded.")

    required_dip_cols = ['Bait', 'string_id']
    if not all(col in dip_id_df.columns for col in required_dip_cols):
         raise KeyError(f"Merged DIP dataframe missing one of required columns: {required_dip_cols}")
    dip_id_df['string_id'] = dip_id_df['string_id'].astype(str)

    dip_proteins_by_bait: Dict[str, List[str]] = {}
    unique_bait_stringid_pairs = dip_id_df[['Bait', 'string_id']].drop_duplicates()
    for bait, group in unique_bait_stringid_pairs.groupby('Bait'):
        dip_proteins_by_bait[bait] = group['string_id'].tolist()
    logging.info(f"Processed DIPs mapped to STRING IDs for {len(dip_proteins_by_bait)} baits.")
    total_dip_ids = sum(len(ids) for ids in dip_proteins_by_bait.values())
    logging.info(f"Total unique (Bait, STRING ID) pairs found for DIPs: {total_dip_ids:,}")

    # Process DEP Data
    if 'identifier' not in dep_df.columns:
        raise KeyError("DEP dataframe ('dep_file') must contain an 'identifier' column holding STRING protein IDs.")
    dep_proteins = dep_df['identifier'].dropna().astype(str).unique().tolist()
    logging.info(f"Extracted {len(dep_proteins):,} unique STRING IDs from DEP data.")

    return dip_proteins_by_bait, dep_proteins


# --- Main Execution Logic ---

def main():
    """
    Main function orchestrating the SIP network construction pipeline,
    modified to process pairs in chunks to manage memory usage.
    """
    # --- 1. Configuration and Setup ---
    try:
        config = load_config("config.json")
        result_dir: Path = config["result_dir"]
        num_processes: int = int(config["num_processes"])
        combined_score_threshold: int = int(config["combined_score"])
        chunk_size: int = int(config["chunk_size"]) # Get chunk size from config
        maxtasksperchild: Optional[int] = config["maxtasksperchild"] # Get maxtasks from config

        result_dir.mkdir(parents=True, exist_ok=True)
        logging.info(f"Results will be saved in: {result_dir}")
        logging.info(f"Using {num_processes} processes for parallel execution.")
        logging.info(f"Processing pairs in chunks of size: {chunk_size:,}")
        if maxtasksperchild:
             logging.info(f"Worker processes will restart after {maxtasksperchild} tasks.")


    except (FileNotFoundError, KeyError, ValueError, TypeError) as e:
         logging.error(f"Configuration error: {e}. Please check config.json. Exiting.")
         sys.exit(1)

    logging.info("Starting Shortest Interaction Path (SIP) network construction.")

    # --- 2. Load Data & Build STRING Graph ---
    try:
        dip_df, dep_df, string_info_df = load_data(
            config["dip_file"], config["dep_file"], config["string_info_file"]
        )
        string_G = build_string_graph(config["string_links_file"], combined_score_threshold)
    except (FileNotFoundError, KeyError, Exception) as e:
        logging.error(f"Failed to load data or build STRING graph: {e}. Exiting.")
        sys.exit(1)

    if string_G.number_of_nodes() == 0 or string_G.number_of_edges() == 0:
        logging.error("The filtered STRING graph is empty. No paths can be calculated. Exiting.")
        sys.exit(1)

    # --- 3. Prepare & Filter Protein Lists ---
    try:
        dip_proteins_by_bait, dep_proteins = get_protein_lists(dip_df, dep_df, string_info_df)
    except KeyError as e:
        logging.error(f"Failed to process protein lists due to missing columns: {e}. Exiting.")
        sys.exit(1)

    nodes_in_graph = set(string_G.nodes())
    logging.info(f"Filtering protein lists against {len(nodes_in_graph):,} nodes present in the STRING graph.")

    dep_proteins_in_graph = [p for p in dep_proteins if p in nodes_in_graph]
    # (Logging for filtered DEPs...)
    if not dep_proteins_in_graph:
         logging.error("No DEP proteins found in the built STRING graph after filtering. Exiting.")
         sys.exit(1)

    dip_proteins_by_bait_in_graph: Dict[str, List[str]] = {}
    # (Filtering and logging for DIPs...)
    total_dips_after_filter = 0
    for bait, dip_list in dip_proteins_by_bait.items():
        filtered_dip_list = [p for p in dip_list if p in nodes_in_graph]
        if filtered_dip_list:
             dip_proteins_by_bait_in_graph[bait] = filtered_dip_list
             total_dips_after_filter += len(filtered_dip_list)
        # else: logging.warning(...)
    if not dip_proteins_by_bait_in_graph:
        logging.error("No DIP proteins found in the built STRING graph after filtering. Exiting.")
        sys.exit(1)
    logging.info(f"{total_dips_after_filter:,} total DIP proteins remain after filtering.")


    # --- 4. Generate Pairs ---
    all_pairs_to_process: List[Tuple[str, str]] = []
    logging.info("Generating (DIP, DEP) pairs...")
    for bait, dip_list_in_graph in dip_proteins_by_bait_in_graph.items():
        bait_pairs = list(itertools.product(dip_list_in_graph, dep_proteins_in_graph))
        bait_pairs_filtered = [(s, t) for s, t in bait_pairs if s != t]
        all_pairs_to_process.extend(bait_pairs_filtered)

    unique_pairs_to_process = list(set(all_pairs_to_process))
    num_total_pairs = len(unique_pairs_to_process)
    logging.info(f"Total unique (DIP, DEP) pairs to process: {num_total_pairs:,}")

    if not unique_pairs_to_process:
        logging.error("No valid pairs found. Exiting.")
        sys.exit(1)

    # --- 5. Calculate Shortest Paths (Chunked Parallel Processing) ---
    final_sip_edges: Set[Tuple[str, str]] = set() # Aggregate edges here directly
    num_chunks = math.ceil(num_total_pairs / chunk_size) # Calculate number of chunks

    logging.info(f"Starting shortest path calculation in {num_chunks} chunks...")

    # Initialize the pool *outside* the loop
    # Add maxtasksperchild here
    pool = mp.Pool(processes=num_processes, maxtasksperchild=maxtasksperchild)
    try:
        for i in range(num_chunks):
            start_index = i * chunk_size
            end_index = min((i + 1) * chunk_size, num_total_pairs)
            chunk_pairs = unique_pairs_to_process[start_index:end_index]

            logging.info(f"Processing chunk {i+1}/{num_chunks} ({len(chunk_pairs):,} pairs)...")

            # Prepare arguments for the current chunk
            chunk_args = [(pair, string_G) for pair in chunk_pairs]

            # Run starmap for the current chunk
            chunk_results = pool.starmap(find_shortest_path_edges, chunk_args)

            # Aggregate results from the chunk immediately
            for edge_list in chunk_results:
                final_sip_edges.update(edge_list)

            logging.info(f"Finished chunk {i+1}/{num_chunks}. Current total unique edges: {len(final_sip_edges):,}")
            # Optional: Add garbage collection hint if memory is extremely tight
            # import gc
            # gc.collect()

    except Exception as e:
        logging.error(f"Multiprocessing failed during chunk processing: {e}")
        # Ensure pool is closed even if an error occurs mid-chunk
        pool.close()
        pool.join()
        sys.exit(1)
    finally:
        # Ensure the pool is always closed and joined
        pool.close()
        pool.join()
        logging.info("Multiprocessing pool closed.")


    # --- 6. Build Final SIP Graph ---
    logging.info("Building final SIP graph from aggregated edges...")
    sip_graph = nx.Graph()
    sip_graph.add_edges_from(list(final_sip_edges)) # Add unique edges collected

    logging.info(f"Constructed final SIP graph: {sip_graph.number_of_nodes():,} nodes, {sip_graph.number_of_edges():,} edges.")

    # (Optional: Add edge attributes back - same as before)

    # --- 7. Save Final Graph ---
    output_file_base = result_dir / f"SIP_network_combined_{combined_score_threshold}"
    try:
        save_obj(sip_graph, output_file_base)
        logging.info(f"Final SIP network graph saved successfully.")
    except Exception as e:
        logging.error(f"Failed to save the final SIP graph to {output_file_base}.pkl: {e}")

    logging.info("SIP network construction process finished.")


# --- Script Entry Point ---
if __name__ == '__main__':
    main()

