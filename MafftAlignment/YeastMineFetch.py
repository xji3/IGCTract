# A file to use YeastMine to get exon intron info
# Xiang Ji
# xji3@ncsu.edu

from intermine.webservice import Service
service = Service("http://yeastmine.yeastgenome.org/yeastmine/service")


if __name__ == '__main__':
    query = service.new_query("ORF")

    query.add_view(
    "organism.name", "organism.taxonId", "CDSs.locations.strand",
    "CDSs.locations.start", "CDSs.locations.end")

    pairs = []
    with open('../All_Pairs.txt', 'r') as f:
        for line in f.readlines():
            pairs.append(line.replace('\n','').split('_'))

    #pairs.remove(['YLR028C', 'YMR120C'])
    genes = [gene for pair in pairs for gene in pair]
    #genes = [pair[0] for pair in pairs]

    genes_feature_list = {}
    ct = 0
    desired_features = ["ORF.CDSs.locations.strand",
    "ORF.CDSs.locations.start", "ORF.CDSs.locations.end"]
    for row in query.rows():
        if row['ORF.secondaryIdentifier'] in genes:
            ct += 1
            print row['ORF.secondaryIdentifier'], ct
            if genes_feature_list.has_key(row['ORF.secondaryIdentifier']):
                genes_feature_list[row['ORF.secondaryIdentifier']].append(dict(row))
            else:
                genes_feature_list[row['ORF.secondaryIdentifier']] = [dict(row)]

    
    gene_to_cds = {gene:sorted([[int(genes_feature_list[gene][i][feat]) for feat in desired_features] for i in range(len(genes_feature_list[gene]))]) for gene in genes}
    for pair in pairs:
        print pair[0], gene_to_cds[pair[0]]
        print pair[1], gene_to_cds[pair[1]]
        print

    with open('./all_pairs_CDS_positions.txt', 'w+') as f:
        for pair in pairs:
            to_write = ['_'.join(pair), str(gene_to_cds[pair[0]][0][0])]
            if gene_to_cds[pair[0]][0][0] == 1:
                start_pos = gene_to_cds[pair[0]][0][1] - 1
                pos_list = [str(abs(item[1] - start_pos)) + '\t' + str(abs(item[2] - start_pos)) for item in gene_to_cds[pair[0]]]
            else:
                start_pos = gene_to_cds[pair[0]][-1][2] + 1
                pos_list = [str(abs(item[2] - start_pos)) + '\t' + str(abs(item[1] - start_pos)) for item in gene_to_cds[pair[0]]]

            if not gene_to_cds[pair[0]][0][0] == 1:
                pos_list.reverse()
            to_write.extend(pos_list)
                
            f.write('\t'.join(to_write) + '\n')
