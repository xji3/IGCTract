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
    with open('../Filtered_pairs.txt', 'r') as f:
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

    
    gene_to_cds = {gene:[[int(genes_feature_list[gene][i][feat]) for feat in desired_features] for i in range(len(genes_feature_list[gene]))] for gene in genes}
    for pair in pairs:
        print pair[0], gene_to_cds[pair[0]]
        print pair[1], gene_to_cds[pair[1]]
        print 
