import hiperwalk as hpw
import networkx as nx
import numpy as np
import pandas as pd
import pickle
import random

seed_value = 777
random.seed(seed_value)
np.random.seed(seed_value) 

long_covid_Graph_df = pd.read_csv('MP_SIP_longcovid_network.csv', names=['node1', 'node2'])
long_covid_Graph_raw = nx.Graph()
long_covid_Graph_raw.add_edges_from(long_covid_Graph_df.values)

id_to_number = dict(zip(long_covid_Graph_raw.nodes(), range(len(long_covid_Graph_raw.nodes()))))
number_to_id = {v: k for k, v in id_to_number.items()}
covid_node = id_to_number['COVID19']
covid_Graph = nx.relabel_nodes(long_covid_Graph_raw, id_to_number)

covid_Graph_adj = nx.adjacency_matrix(covid_Graph)
covid_Graph_hpw = hpw.Graph(covid_Graph_adj)

print('Hyperwalk Graph Statistics')
print('number_of_vertices', covid_Graph_hpw.number_of_vertices())
print('number_of_edges', covid_Graph_hpw.number_of_edges())
print('number_of_loops', covid_Graph_hpw.number_of_loops())
print('vertex_number', covid_Graph_hpw.vertex_number(covid_node))

pagerank_result = nx.pagerank(covid_Graph, personalization={covid_node: 1}, max_iter=10000, alpha=0.85)
pagerank_result_sr = pd.Series(pagerank_result)
pagerank_result_sr.rename(index= number_to_id, inplace=True)

pickle.dump(pagerank_result_sr, open('/data/jp2079/QWP_Data/result/COVID_RWR_result.pkl', 'wb'))

coined = hpw.Coined(graph = covid_Graph_hpw)
coined_uniform_state = coined.uniform_state([covid_node])
coined_final_state = coined.simulate(range=(10000),state=coined_uniform_state)

coined_result = coined.probability_distribution(coined_final_state)
coined_result_df = pd.DataFrame(coined_result).T
coined_result_df.rename(index= number_to_id, inplace=True)

pickle.dump(coined_result_df, open('COVID_coined_result.pkl, 'wb'))