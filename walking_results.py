import numpy as np
import pandas as pd
import pickle
import networkx as nx

from scipy.spatial.distance import euclidean

rwr_result_df = pickle.load(open('COVID_RWR_result.pkl', 'rb'))
coined_result_df = pickle.load(open('COVID_coined_result.pkl', 'rb'))

covid_Graph_df = pd.read_csv('MP_SIP_longcovid_network.csv', names=['node1', 'node2'])
covid_Graph = nx.Graph()
covid_Graph.add_edges_from(covid_Graph_df.values)

protein_df = pd.read_table('9606.protein.info.v12.0.txt')
protein_df['string_id']= protein_df['#string_protein_id'].str.split('.').str[1]
protein_df= protein_df.set_index('string_id')

def AveragedProbabilityDistribution(df):
    summ = df.loc[:,0]
    average_dict = {0:summ}
    for i in range(1, df.shape[1]):
        summ = summ + df.loc[:,i]
        average = summ / (i+1)
        average_dict[i] = average
        
    return pd.DataFrame(average_dict)

def MixingTime(target):
    for i,v in enumerate(target):
        if v < 0.00001:
            return i
            
coined_averaged_df = AveragedProbabilityDistribution(coined_result_df)
coined_averaged_sim = [euclidean(coined_averaged_df[i], coined_averaged_df[i+1]) for i in range(9999)]
coined_mixing_time = MixingTime(coined_averaged_sim)

def FinalDataGeneration():
    shortest_result = list()
    for i in covid_Graph.nodes():
        length = nx.shortest_path_length(covid_Graph,'COVID19',i)
        shortest_result.append({'target':i, 'distance':length, 'rwr':rwr_result_df.loc[i],
                                'coined': coined_averaged_df.loc[i,coined_mixing_time]})
    shortest_df = pd.DataFrame.from_dict(shortest_result)
    shortest_df.set_index('target',inplace=True)

    merged_df = pd.concat([coined_averaged_df[coined_mixing_time], rwr_result_df],axis=1)
    
    merged_df.columns=['coined','rwr']
    merged_df = merged_df.join(protein_df['preferred_name'])
    merged_df = merged_df.join(shortest_df['distance'])

    return merged_df

result_df = FinalDataGeneration()
result_df.to_csv('result.csv')
