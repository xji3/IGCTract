from IGCexpansion.CodonGeneconv import ReCodonGeneconv
import argparse

def main(args):
    paralog = [args.paralog1, args.paralog2]
    Force = None
    alignment_file = './MafftAlignment/' + '_'.join(paralog) + '/' + '_'.join(paralog) + '_input.fasta'
    newicktree = './YeastTree.newick'
    if args.force:
        if args.model == 'MG94':
            Force = {5:0.0}
        elif args.model == 'HKY':
            Force = {4:0.0}
    else:
        Force = None
    test = ReCodonGeneconv( newicktree, alignment_file, paralog, Model = args.model, Force = Force, clock = args.clock)
    test.get_mle(True, True, 0, 'BFGS')
    test.get_individual_summary(summary_path = './summary/')
    #test.get_SitewisePosteriorSummary(summary_path = './Summary/')
    if Force == None:
        test.get_sitewise_loglikelihood_summary('./summary/' + '_'.join(paralog) + '_MG94_nonclock_sw_lnL.txt')
    else:
        test.get_sitewise_loglikelihood_summary('./summary/Force_' + '_'.join(paralog) + '_MG94_nonclock_sw_lnL.txt')


if __name__ == '__main__':
##    parser = argparse.ArgumentParser()
##    parser.add_argument('--model', required = True, help = 'Substitution Model')
##    parser.add_argument('--paralog1', required = True, help = 'Name of the 1st paralog')
##    parser.add_argument('--paralog2', required = True, help = 'Name of the 2nd paralog')
##    parser.add_argument('--force', dest = 'force', action = 'store_true', help = 'Tau parameter control')
##    parser.add_argument('--no-force', dest = 'force', action = 'store_false', help = 'Tau parameter control')
##    parser.add_argument('--clock', dest = 'clock', action = 'store_true', help = 'clock control')
##    parser.add_argument('--no-clock', dest = 'clock', action = 'store_false', help = 'clock control')
##    
##    main(parser.parse_args())

    pairs = []
    all_pairs = './Filtered_pairs.txt'
    with open(all_pairs, 'r') as f:
        for line in f.readlines():
            pairs.append(line.replace('\n','').split('_'))

    for pair in pairs[:1]:
        paralog = pair
        Force = None
        #Force = {5:0.0}
        alignment_file = './MafftAlignment/' + '_'.join(paralog) + '/' + '_'.join(paralog) + '_input.fasta'
        newicktree = './YeastTree.newick'
        test = ReCodonGeneconv( newicktree, alignment_file, paralog, Model = 'MG94', Force = Force, clock = None)
        test.get_mle(True, True, 0, 'BFGS')
        test.get_sitewise_loglikelihood_summary('./summary/' + '_'.join(paralog) + '_MG94_nonclock_sw_lnL.txt')

        Force = {5:0.0}
        test_force = ReCodonGeneconv( newicktree, alignment_file, paralog, Model = 'MG94', Force = Force, clock = None)
        test_force.update_by_x(test.x)
        test_force._loglikelihood2()
        test_force.get_sitewise_loglikelihood_summary('./summary/Force_' + '_'.join(paralog) + '_MG94_nonclock_sw_lnL.txt')

