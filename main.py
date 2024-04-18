import sys
import argparse
import time
from Functions.construct_mutation_matrix import *
from Functions.construct_bipartite_graph import *
from Functions.LBSA import *
from Functions.results import *


def main():
    #Here we read the necessary parameters required for running the algorithm
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, \
                                    usage='\n\npython main.py <-n network file> <-m dataset name> [-c cooling factor] [-i number of iterations]\n', \
                                    description='', epilog='Implementation of MYTHESIS.\n')

    parser.add_argument('-n', '--network', required=True, type=str, help='mutations file.', metavar='', dest="network")
    parser.add_argument('-d', '--dataset', required=True, type=str, help='dataset name.', metavar='', dest="dataset_name")
    parser.add_argument('-c', '--cooling', default=(1 - (10 ** -2)), type=int, help='Cooling factor', metavar='', dest="cooling")
    parser.add_argument('-i', '--iterations', default=(10 ** 5), type=int, help='Number of iterations', metavar='', dest="iterations")

    args = parser.parse_args() #all input arguments are parsed into args variable

    lbsa_starttime = time.time()
    #get input data
    maf_file_name = args.network
    dataset_name = args.dataset_name
    cooling_fact = args.cooling
    iterations = args.iterations

    #read data
    print "START RUN...\n Reading input data..."
    mutation_matrix = get_mutation_matrix_from_maf("Data/TCGA/" + maf_file_name)
    influence_graph = nx.read_edgelist("Data/Reactome_FIsInGene_2021.txt", delimiter='\t')
    gene_expression_matrix = pd.read_csv('Data/GBMPatientOutlierMatrix.csv', index_col=0)
    print "Finished Reading input."

    #create bipartite graph
    print "Creating Bipartite Graph. Start..."
    bipartite_graph, green_nodes = create_bipartite_graph(mutation_matrix, influence_graph, gene_expression_matrix)
    print "Bipartite Graph Created Successfully."

    #implement List Based Simulated Annealing
    print "Apply List Based Simulated Annealing"
    lbsa_drivers_list, drivers_order_list2 = List_Based_Simulated_Annaeling(dataset_name, green_nodes, bipartite_graph, cooling_fact, iterations)
    print "Algorithm successful..."

    #write output
    print "Writing Data to output..."
    report(lbsa_drivers_list, drivers_order_list2, dataset_name)
    print "Output file -- ranked_driver_genes_" + dataset_name + ".txt -- created.\n Check Output Folder"
    print "Cheers!!!"

    lbsa_endtime = time.time()
    total_algorithm_runtime = lbsa_endtime - lbsa_starttime
    total_algorithm_runtime_mins = total_algorithm_runtime / 60
    print("---Total program runtime ----%s minutes ----" % (total_algorithm_runtime_mins))

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        print "Error Encountered. EXITING...."
        sys.exit(0)

