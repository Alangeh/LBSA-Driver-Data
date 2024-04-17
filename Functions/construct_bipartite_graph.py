import networkx as nx


def patient_mutation_exist(gene, patient, mutation_matrix):
  return (mutation_matrix.at[patient,gene] == 1)

def verify_interaction(gi, gj, influence_graph):
  return influence_graph.has_edge(gi, gj)

def verify_gene_expression(gene, patient, gene_expression_matrix):
  return gene_expression_matrix.at[patient,gene]


def create_bipartite_graph(mutation_matrix, influence_graph, gene_expression_matrix):
    bipartite_graph = nx.Graph()
    patients = mutation_matrix.index.to_list()
    genes = mutation_matrix.columns.to_list();

    left_partition = genes
    bipartite_graph.add_nodes_from(left_partition, bipartite=0)

    right_partition = []
    for p in patients:
        for gene in genes:
            right_partition.append(gene + "_" + p)
    bipartite_graph.add_nodes_from(right_partition, bipartite=1)

    green_nodes = set()
    for gi in left_partition:
        for pk in patients:
            if patient_mutation_exist(gi, pk, mutation_matrix):
                green_nodes.add(gi)
                for gj in genes:
                    if verify_interaction(gi, gj, influence_graph):
                        #if (verify_gene_expression(gj, pk, gene_expression_matrix)):
                        label_node_right_partition = gj + "_" + pk
                        bipartite_graph.add_edge(gi, label_node_right_partition)



    return bipartite_graph, green_nodes