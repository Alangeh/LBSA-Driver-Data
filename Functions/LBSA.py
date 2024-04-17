import random
import math
import numpy as np
from results import *
Driver_Count = 100

def calc_coverage(neighbors_lists):
  all_neighbors_set = set()
  for neighbors_list in neighbors_lists:
    all_neighbors_set.update(neighbors_list)
  return len(all_neighbors_set)

def get_parameters(t, c, i):
    temp = t
    cool_f = c
    iterations = i
    return temp, cool_f, iterations

def compute_p(drivers=None, random_res=None):

    pvalues_dict = {key:1 for key in drivers}

    tmp_dict = {key:0 for key in drivers}
    total_random_drivers = 0
    for i in random_res: # compare results of randomly generated drivers with drivers from original data
        total_random_drivers += len(random_res[i])
        for driver in drivers:
            for random_driver in random_res[i]:
                if random_res[i][random_driver] > drivers[driver]:
                    tmp_dict[driver] += 1

    for driver in tmp_dict: #calculate the p - values for each driver in generated drivers dictionary
        pvalues_dict[driver] = float(tmp_dict[driver])/total_random_drivers

    return pvalues_dict


def Generate_Temperature_List():
    f = 100
    x = 99.0
    cool_factor = 10 ** 6
    temperature_list = []
    temperature_list1 = []
    max_list_length = 10
    p0 = 0.1

    for i in range(0, 1000):
        y = float(random.randint(2, f))
        f_x = 1 / x
        f_y = 1 / y
        if (f_y < f_x):
            x = y
        x = x-1
        t = -1 * (abs(f_y - f_x) / np.log(p0))
        t = int(t * cool_factor)
        if (t < 1000):
            temperature_list.append(t)

    temp_list = []
    temperature_list1 = list(set(temperature_list))
    temperature_list1.sort(reverse=True)
    for j in range (max_list_length):
        temp_list.append(temperature_list1[j])

    return temp_list

def List_Based_Simulated_Annaeling(dataset_name, green_nodes, bipartite_graph, cool=1 - (10 ** -2), iterations=10 ** 5):
    neighbors_list = []
    temperature_list = []
    temperature_list = Generate_Temperature_List()
    temp = max(temperature_list)
    new_temperature = 0
    for g in green_nodes:
        neighbors_list.append(list(bipartite_graph.neighbors(g)))
    neighbors_dictionary = dict(zip(green_nodes, neighbors_list))
    print 'Define initial parameters for list based simulated annealing'
    initial_temperature, cooling_factor, iterations = get_parameters(temp, cool, iterations)

    left_partition_degree = list(green_nodes)

    best_set = []
    random_number_list = random.sample(range(0, len(left_partition_degree) - 1), Driver_Count)
    for i in range(0, Driver_Count):
        g = left_partition_degree[random_number_list[i]]
        best_set.append(g)

    neighbors_lists_best_set = []

    print 'Create best set of genes...'
    for g in best_set:
        neighbors_lists_best_set.append(list(bipartite_graph.neighbors(g)))
    best_coverage = calc_coverage(neighbors_lists_best_set)
    current_temperature = initial_temperature
    print 'Start Iterations...'
    for i in range(0, iterations):
        current_set = list(best_set)  # alternative to .copy()
        # current_set = copy.deepcopy(best_set)
        neighbors_lists_current_set = list(neighbors_lists_best_set)  # alternative to .copy()

        g_j = left_partition_degree[random.randint(0, len(left_partition_degree) - 1)]
        if (not (g_j in best_set)):
            index = random.randint(0, len(current_set) - 1)
            del current_set[index]
            del neighbors_lists_current_set[index]
            current_set.append(g_j)
            neighbors_lists_current_set.append(neighbors_dictionary.get(g_j))

        if (current_set != best_set):
            current_coverage = calc_coverage(neighbors_lists_current_set)
            if (current_coverage > best_coverage):
                best_set = current_set
                best_coverage = current_coverage
                neighbors_lists_best_set = neighbors_lists_current_set
                new_temperature = current_temperature
            else:
                diff = abs(current_coverage - best_coverage)
                if (random.random() < math.exp(-diff / current_temperature)):
                    best_set = current_set
                    best_coverage = current_coverage
                    neighbors_lists_best_set = neighbors_lists_current_set
                    temperature_list = Generate_Temperature_List()
                    new_temperature = max(temperature_list)
                else:
                    new_temperature = current_temperature
        else:
            new_temperature = current_temperature
        current_temperature = new_temperature * cooling_factor
    print 'iterations end...'
    print 'Best set created...'
    print 'Best set count --> ' + str(len(best_set))
    report_drivers(best_set, dataset_name)
    drivers_list1 = []
    print 'Creating ordered list for drivers...'
    for g in best_set:
        drivers_list1.append((g, bipartite_graph.degree[g]))
    drivers_order_list1 = sorted(drivers_list1, key=lambda x: x[1], reverse=True)
    drivers_order_list1 = list(map(lambda x: x[0], drivers_order_list1))

    drivers_order_list2 = []
    print 'Remove neighbours, and calculate p_values'
    for g in drivers_order_list1:
        if (bipartite_graph.degree[g] >= 1 ):
            p_value = 10 / float(bipartite_graph.degree[g])
            drivers_order_list2.append((g, bipartite_graph.degree[g], p_value))

        g_neighbors = list(bipartite_graph.neighbors(g))
        for neighbor in g_neighbors:
            bipartite_graph.remove_node(neighbor)
    drivers_order_list2 = sorted(drivers_order_list2, key=lambda x: x[1], reverse=True)
    lbsa_drivers_list = list(map(lambda x: x[0], drivers_order_list2))
    print 'Drivers list successfully generated!!!'

    return lbsa_drivers_list, drivers_order_list2