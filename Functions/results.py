#write output to file, out put file consist of ranked driver genes.

def report(lbsa_drivers_list, drivers_order_list2, dataset_val):

    output_prefix = 'ranked_driver_genes_' + dataset_val
    output = open('Output/%s.txt' % output_prefix, 'w')
    output.write('\t'.join(['rank', 'degree', 'gene', 'p-value']) + '\n')
    print ('rank\tdegree\tgene\tp-value\n')
    for ix, item in enumerate(drivers_order_list2):
        gene = item[0]
        degree = item[1]
        p_value = item[2]
        output.write('\t'.join([str(ix+1), str(degree), str(gene), str(p_value)]) + '\n')
        print(str(ix+1) + '\t' + str(degree) + '\t' + str(gene) + '\t' + str(p_value) + '\n')
    output.close()

def report_drivers(best_set, dataset_name):
    output_prefix = 'drivers_list_' + dataset_name
    output = open('Output/%s.txt' % output_prefix, 'w')
    output.write('\t'.join(['No.', 'DRIVER GENES']) + '\n')
    for ix, item in enumerate(best_set):
        output.write('\t'.join([str(ix+1), str(item)]) + '\n')
    output.close()
#test
#list = ['a']
#list1 = [('a', 1, 0.1),('b', 3, 0.01),('c',2,0.2)]
#report(list,list1,'prad')