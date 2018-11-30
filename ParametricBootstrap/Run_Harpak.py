# This is a separate file to run all 20 Harpak data sets together
# This file implements a naiive parallel computing method without careful testing
# Xiang Ji
# xji3@ncsu.edu


# OK, naiive parallel computing upgrade
# http://sebastianraschka.com/Articles/2014_multiprocessing.html

from IGCexpansion.PSJSGeneconv import PSJSGeneconv
from IGCexpansion.JSGeneconv import JSGeneconv
import numpy as np
import scipy, os, argparse
import multiprocessing as mp
from functools import partial
from collections import namedtuple
from copy import deepcopy

# pickle dependent
# https://stackoverflow.com/questions/1816958/cant-pickle-type-instancemethod-when-using-multiprocessing-pool-map/7309686#7309686
##
##from copy_reg import pickle
##from types import MethodType
##
##def _pickle_method(method):
##    func_name = method.im_func.__name__
##    obj = method.im_self
##    cls = method.im_class
##    return _unpickle_method, (func_name, obj, cls)
##
##def _unpickle_method(func_name, obj, cls):
##    for cls in cls.mro():
##        try:
##            func = cls.__dict__[func_name]
##        except KeyError:
##            pass
##        else:
##            break
##    return func.__get__(obj, cls)

class Run_PSJS_Harpak_all:
    auto_save_step = 2
    def __init__(self, alignment_file_list, gene_to_orlg_file,
                 seq_index_file_list, 
                 tree_newick, DupLosList,
                 x_js, pm_model, IGC_pm, rate_variation,
                 node_to_pos, terminal_node_list,
                 save_file_list, log_file_list,
                 save_file,
                 multiprocess_combined_list = None,
                 num_processes = 1,                                # default number of parallel processes
                 cdna = False, allow_same_codon = False,
                 force = None, nsites = None, space_list = None):
        # first check if the length of all input lists agree
        assert(len(set([len(alignment_file_list), len(seq_index_file_list), \
                        len(save_file_list), len(log_file_list)])) == 1)
        self.psjsgeneconv_list = [PSJSGeneconv(alignment_file_list[i], gene_to_orlg_file, seq_index_file_list[i],
                                               False, False, tree_newick, DupLosList, x_js, pm_model, IGC_pm,
                                               False, node_to_pos, terminal_node_list, save_file_list[i], log_file_list[i])\
                                  for i in range(len(alignment_file_list))]
        self.save_file = save_file
        self.x = None
        self.num_processes = num_processes
        self.multiprocess_combined_list = multiprocess_combined_list

        if os.path.isfile(self.save_file):
            self.initialize_by_save()
            print ('Loaded parameters from ' + self.save_file)
        else:
            self.x   = np.concatenate((x_js, np.array([-0.69682578, -1.45309823, -3.12708232, -4.15295349, -4.93252404,
                                                       -3.97566549, -3.90817818, -4.71788783, -4.49238764, -4.54912794])))
        self.auto_save = 0

    def unpack_x(self, x):
        self.x = x
        for psjsgeneconv in self.psjsgeneconv_list:
            psjsgeneconv.unpack_x(x)

    def _process_objective_and_gradient(self, num_psjsgeneconv_list, display, x, output):
        for num_psjsgeneconv in num_psjsgeneconv_list:
            result = self.psjsgeneconv_list[num_psjsgeneconv].objective_and_gradient(display, x)
            output.put(result)


    def objective_and_gradient(self,display, x, parallel = True):
        self.unpack_x(deepcopy(x))

        f = 0.0
        g = 0.0
        if parallel:            
            # Define an output queue
            output = mp.Queue()
            
            # Setup a list of processes that we want to run
            if self.multiprocess_combined_list is None:
                processes = [mp.Process(target=self._process_objective_and_gradient, args=([i], False, x, output)) \
                             for i in range(len(self.psjsgeneconv_list))]
            else:
                processes = [mp.Process(target=self._process_objective_and_gradient, args=(i, False, x, output)) \
                             for i in self.multiprocess_combined_list]
            
            # Run processes
            for p in processes:
                p.start()

            # Exit the completed processes
            for p in processes:
                p.join()

            # Get process results from the output queue
            results = [output.get() for p in processes]

    ##        pool = mp.Pool(processes = self.num_processes)
    ##        results = [pool.apply(psjsgeneconv.objective_and_gradient, args = (display, x))\
    ##                   for psjsgeneconv in self.psjsgeneconv_list]

            for result in results:
                f += result[0]
                g += result[1]
        else:
            for psjsgeneconv in self.psjsgeneconv_list:
                f_inc, g_inc = psjsgeneconv.objective_and_gradient(False, x)
                f += f_inc
                g += g_inc

        if display:
            print('log likelihhood = ' ,f)
            print('exp x = ',  np.exp(self.x))
            print('Gradient = ', g)
                
        self.perform_auto_save()
        return f, g

    def save_x(self):
        for psjsgeneconv in self.psjsgeneconv_list:
            psjsgeneconv.save_x()
        save = self.x
        np.savetxt(self.save_file, save.T)

    def initialize_by_save(self):
        self.x = np.loadtxt(open(self.save_file, 'r'))
        self.unpack_x(self.x)
            
    def objective_wo_gradient(self, display, x):
        self.unpack_x(x)
        f = 0.0
        for psjsgeneconv in self.psjsgeneconv_list:
            f_incr = self.objective_wo_gradient(display = display, x = x)
            f += f_incr
        return f

    def perform_auto_save(self):
        # Now save parameter values
        self.auto_save += 1
        if self.auto_save == Run_PSJS_Harpak_all.auto_save_step:
            self.save_x()
            self.auto_save = 0

    def get_gradient_hessian(self, x, gradient_file_list, hessian_file_list):
        for i, psjsgeneconv in enumerate(self.psjsgeneconv_list):

            psjsgeneconv.get_gradient_hessian(x, gradient_file_list[i], hessian_file_list[i])

            if os.path.isfile(gradient_file_list[i]) and os.path.isfile(hessian_file_list[i]):
                continue
            else:
                psjsgeneconv.get_gradient_hessian(x, gradient_file_list[i], hessian_file_list[i])

    
    def get_mle(self, display = True, derivative = True, stringent_level = 'low'):
        self.unpack_x(self.x)  # do one more update first
        if derivative:
            f = partial(self.objective_and_gradient, display)
        else:
            f = partial(self.objective_wo_gradient, display)

        if stringent_level == 'low':
            factr = 1e12
        elif stringent_level == 'moderate':
            factr = 1e7
        elif stringent_level == 'high':
            factr = 10.0
        else:
            exit('Check stringent_level in get_mle() function!')

        guess_x = self.x
        
        # only implemented 'One rate' case of IGC rate parameterization
        # This is a lazy constraint that should be removed later
        assert(len(self.psjsgeneconv_list[0].psjsmodel.x_js) == len(self.psjsgeneconv_list[0].psjsmodel.x_pm) + 2)
        
        bnds = [(None, -0.001)] * 3
        bnds.extend([(None, None)] * (len(self.psjsgeneconv_list[0].psjsmodel.x_pm) - 3 + 1))
        bnds.extend([(None, 0.0)])
        bnds.extend([(None, None)] * (len(self.x) - len(self.psjsgeneconv_list[0].psjsmodel.x_js)))

##        if self.root_by_dup:
##            bnds  = [(None, None)] * len(self.tree.edge_list)
##        else:
##            bnds  = [(None, None)] * (len(self.tree.edge_list) - 1)



        if derivative:
            result = scipy.optimize.minimize(f, guess_x, jac = True, method = 'L-BFGS-B', bounds = bnds, options={'ftol':factr*1e-08})
        else:
            result = scipy.optimize.minimize(f, guess_x, jac = False, method = 'L-BFGS-B', bounds = bnds, options={'ftol':factr*1e-08})

        self.save_x()
        print(result)
        return result

    def objective_tract_p(self, display, log_tract_p):
        f = 0.0
        self.update_only_log_tract_p(log_tract_p)
        for psjsgeneconv in self.psjsgeneconv_list:
            f_inc = psjsgeneconv.objective_tract_p(True, log_tract_p)
            f += f_inc
        if display:
            print('log likelihhood = ' ,f)
            print('exp x = ',  np.exp(self.x))
        self.perform_auto_save()
        return f

    def update_only_log_tract_p(self, log_tract_p):
        self.x[4] = self.x[4] - self.x[5] + log_tract_p
        self.x[5] = log_tract_p        

    def optimize_x_IGC(self, guess, display = True, dimension = 1, method = 'L-BFGS-B'):

        f = partial(self.objective_tract_p, display)
        guess_x = guess
        if method == 'L-BFGS-B':
            bnds = [(None, 0.0)]
        elif method == 'BasinHopping':
            bnds = [(-20.0, 0.0)]
        elif method == 'DifferentialEvolution':
            bnds = [(-20.0, 0.0)]
        else:
            sys.exit('Optimization method is not implemented!')

        if method == 'L-BFGS-B':
            result = scipy.optimize.minimize(f, guess_x, jac = False, method = method, bounds = bnds)
        else:
            sys.exit('Not implemented yet!')
            
        self.save_x()
        print (result)
        return result

class Run_JS_Harpak_all:
    auto_save_step = 5
    def __init__(self, alignment_file_list, gene_to_orlg_file,
                 seq_index_file_list, 
                 tree_newick, DupLosList,
                 x_js, pm_model, IGC_pm, rate_variation,
                 node_to_pos, terminal_node_list,
                 save_file_list, log_file_list,
                 save_file,
                 num_processes = 1,                                # default number of parallel processes
                 cdna = False, allow_same_codon = False,
                 force = None, nsites = None, space_list = None):
        # first check if the length of all input lists agree
        assert(len(set([len(alignment_file_list), len(seq_index_file_list), \
                        len(save_file_list), len(log_file_list)])) == 1)
        self.jsgeneconv_list = [JSGeneconv(alignment_file_list[i], gene_to_orlg_file, False,
                                           tree_newick, DupLosList, x_js, pm_model, IGC_pm, False,
                                           node_to_pos, terminal_node_list,
                                            save_file_list[i], force = force)\
                                  for i in range(len(alignment_file_list))]
        self.save_file = save_file
        self.x = None
        self.num_processes = num_processes
        self.auto_save = 0

        if os.path.isfile(self.save_file):
            self.initialize_by_save()
            print ('Loaded parameters from ' + self.save_file)
        else:
            self.x   = np.array([-1.80846702, -0.24314432, -0.44540066,  0.44970288, 2.99036552, \
                                 0.49424267, -4.57951986, -2.96671004, -4.65348427,\
                                 -4.41268635, -4.41431052, -4.83975226, -4.84754018, -3.59631986,\
                                 -4.82976336])
            self.auto_save = 0

    def unpack_x(self, x):
        self.x = x
        for jsgeneconv in self.jsgeneconv_list:
            jsgeneconv.unpack_x(x)

    def _process_objective_and_gradient(self, num_jsgeneconv, display, x, output):
        result = self.jsgeneconv_list[num_jsgeneconv].objective_and_gradient(display, x)
        output.put(result)

    def objective_and_gradient(self,display, x, parallel = False):
        self.unpack_x(x)

        f = 0.0
        g = 0.0
        if parallel:
            # Define an output queue
            output = mp.Queue()
            
            # Setup a list of processes that we want to run
            processes = [mp.Process(target=self._process_objective_and_gradient, args=(i, False, x, output)) \
                         for i in range(len(self.jsgeneconv_list))]
            
            # Run processes
            for p in processes:
                p.start()

            # Exit the completed processes
            for p in processes:
                p.join()

            # Get process results from the output queue
            results = [output.get() for p in processes]

    ##        pool = mp.Pool(processes = self.num_processes)
    ##        results = [pool.apply(jsgeneconv.objective_and_gradient, args = (display, x))\
    ##                   for jsgeneconv in self.jsgeneconv_list]

            for result in results:
                f += result[0]
                g += result[1]
        else:
            for jsgeneconv in self.jsgeneconv_list:
                f_inc, g_inc = jsgeneconv.objective_and_gradient(False, x)
                f += f_inc
                g += g_inc

        if display:
            print('log likelihhood = ' ,f)
            print('exp x = ',  np.exp(self.x) )
            print('Gradient = ', g)

        self.perform_auto_save()
            
        return f, g

    def perform_auto_save(self):
        # Now save parameter values
        self.auto_save += 1
        if self.auto_save == Run_JS_Harpak_all.auto_save_step:
            self.save_x()
            self.auto_save = 0

    def save_x(self):
        for jsgeneconv in self.jsgeneconv_list:
            jsgeneconv.save_x()
        save = self.x
        np.savetxt(self.save_file, save.T)

    def initialize_by_save(self):
        self.x = np.loadtxt(self.save_file)
        self.unpack_x(self.x)
            
    def objective_wo_gradient(self, display, x):
        self.unpack_x(x)
        f = 0.0
        for jsgeneconv in self.jsgeneconv_list:
            f_incr = self.objective_wo_gradient(display = display, x = x)
            f += f_incr
        return f


    
    def get_mle(self, display = True, derivative = True, stringent_level = 'low'):
        self.unpack_x(self.x)  # do one more update first
        if derivative:
            f = partial(self.objective_and_gradient, display)
        else:
            f = partial(self.objective_wo_gradient, display)

        if stringent_level == 'low':
            factr = 1e12
        elif stringent_level == 'moderate':
            factr = 1e7
        elif stringent_level == 'high':
            factr = 10.0
        else:
            exit('Check stringent_level in get_mle() function!')

        guess_x = self.x
        
        # only implemented 'One rate' case of IGC rate parameterization
        # This is a lazy constraint that should be removed later
        assert(len(self.jsgeneconv_list[0].jsmodel.x_js) == len(self.jsgeneconv_list[0].jsmodel.x_pm) + 1)
        
        bnds = [(None, -0.001)] * 3
        bnds.extend([(None, None)] * (len(self.jsgeneconv_list[0].jsmodel.x_pm) - 3 + 1))
        bnds.extend([(None, None)] * (len(self.x) - len(self.jsgeneconv_list[0].jsmodel.x_js)))

##        if self.root_by_dup:
##            bnds  = [(None, None)] * len(self.tree.edge_list)
##        else:
##            bnds  = [(None, None)] * (len(self.tree.edge_list) - 1)



        if derivative:
            result = scipy.optimize.minimize(f, guess_x, jac = True, method = 'L-BFGS-B', bounds = bnds, options={'ftol':factr*1e-08})
        else:
            result = scipy.optimize.minimize(f, guess_x, jac = False, method = 'L-BFGS-B', bounds = bnds, options={'ftol':factr*1e-08})

        self.save_x()
        print(result)
        return result




def main(args):
    seq_file_list = np.loadtxt('../Harpak/missing_0_species_list.txt', dtype = str)
    pairs = [file_name.split(".")[:2] for file_name in seq_file_list]
    sim_num = args.sim_num
    case_folder = args.case

    alignment_file_list = ['./bootstrapReplicates/' + case_folder + '/' + '_'.join(paralog) \
                           + '/sim_' + str(sim_num) + '/' + '_'.join(paralog) + '_sim_' + str(sim_num) + '.fasta' for paralog in pairs]
    seq_index_file_list = ['../Harpak/prepared_input/' + seq_file.replace('.pos.seq.formatted', '').replace('.', '_') \
                           + '_seq_index.txt' for seq_file in seq_file_list]


    gene_to_orlg_file = '../Harpak/GeneToOrlg.txt'
    save_file = './save/' + case_folder + '/JS_Grand_save_sim_' + str(sim_num) + '.txt'

    tree_newick = '../Harpak/HarpakTree.newick'
    DupLosList = '../Harpak/HarpakDupLost.txt'
    terminal_node_list = ['Human', 'Chimp', 'Goril', 'Orang', 'Macaq', 'Mouse']
    node_to_pos = {'D1':0}

    ###### Now get HKY+PSJS-IGC estimates
    IGC_pm = 'One rate'
    pm_model = 'HKY'
    for paralog in pairs:
        log_paralog_folder = './log/' + case_folder + '/' + '_'.join(paralog)
        log_sim_folder = log_paralog_folder + '/sim_' + str(sim_num)
        save_paralog_folder = './save/' + case_folder + '/' + '_'.join(paralog)
        save_sim_folder = save_paralog_folder + '/sim_' + str(sim_num)
        summary_paralog_folder = './summary/' + case_folder + '/' + '_'.join(paralog)
        summary_sim_folder = summary_paralog_folder + '/sim_' + str(sim_num)

        if not os.path.isdir(log_paralog_folder):
            os.mkdir(log_paralog_folder)
        if not os.path.isdir(log_sim_folder):
            os.mkdir(log_sim_folder)
        if not os.path.isdir(save_paralog_folder):
            os.mkdir(save_paralog_folder)
        if not os.path.isdir(save_sim_folder):
            os.mkdir(save_sim_folder)
        if not os.path.isdir(summary_paralog_folder):
            os.mkdir(summary_paralog_folder)
        if not os.path.isdir(summary_sim_folder):
            os.mkdir(summary_sim_folder)

    save_file_list = ['./save/' + case_folder + '/' + '_'.join(paralog) + '/sim_' + str(sim_num) + '/JS_' + '_'.join(paralog) \
                           + '_guess_' + str(args.guess) + '_save.txt' for paralog in pairs]
    log_file_list = ['./log/' + case_folder + '/' + '_'.join(paralog) + '/sim_' + str(sim_num) + '/JS_' + '_'.join(paralog) \
                           + '_guess_' + str(args.guess) + '_log.txt' for paralog in pairs]
    summary_file_list = ['./summary/' + case_folder + '/' + '_'.join(paralog) + '/sim_' + str(sim_num) + '/JS_' + '_'.join(paralog) \
                           + '_guess_' + str(args.guess) + '_summary.txt' for paralog in pairs]

    rate_variation = False


    x_js = np.log([0.3,0.4,0.5,4.0, 3.0])
      
    test = Run_JS_Harpak_all(alignment_file_list, gene_to_orlg_file,
                 seq_index_file_list, 
                 tree_newick, DupLosList,
                 x_js, pm_model, IGC_pm, rate_variation,
                 node_to_pos, terminal_node_list,
                 save_file_list, log_file_list, save_file, 10)
    test.get_mle(stringent_level = 'high')

    save_file = './save/' + case_folder + '/PSJS_Grand_save_sim_' + str(sim_num) + '.txt'

    save_file_list = ['./save/' + case_folder + '/' + '_'.join(paralog) + '/sim_' + str(sim_num) + '/PSJS_' + '_'.join(paralog) \
                           + '_guess_' + str(args.guess) + '_save.txt' for paralog in pairs]
    log_file_list = ['./log/' + case_folder + '/' + '_'.join(paralog) + '/sim_' + str(sim_num) + '/PSJS_' + '_'.join(paralog) \
                           + '_guess_' + str(args.guess) + '_log.txt' for paralog in pairs]
    summary_file_list = ['./summary/' + case_folder + '/' + '_'.join(paralog) + '/sim_' + str(sim_num) + '/PSJS_' + '_'.join(paralog) \
                           + '_guess_' + str(args.guess) + '_summary.txt' for paralog in pairs]
                           

    initial_tract_length_list = np.log([30.0, 200.0, 500.0])
    guess_lnp = -initial_tract_length_list[args.guess - 1]

    rate_variation = False
    true_x = np.loadtxt('./trueValues/PSJS_HKY_Harpak_One_rate_nonclock_save.txt')
    x_js = np.concatenate((true_x[:4], [true_x[4] - true_x[5] + guess_lnp, guess_lnp]))
    guess_x = np.concatenate((x_js, true_x[6:])) 


    # multiprocess_combined_list = [[0, 9, 12], [7, 5], [6, 1], [18, 8],
    #                               [2], [3], [4], [10], [11],
    #                               [13], [14], [15], [16], [17],
    #                               [19]]

    multiprocess_combined_list = None
      
    test = Run_PSJS_Harpak_all(alignment_file_list, gene_to_orlg_file,
                               seq_index_file_list,
                               tree_newick, DupLosList,
                               x_js, pm_model, IGC_pm, rate_variation,
                               node_to_pos, terminal_node_list,
                               save_file_list, log_file_list, save_file,
                               multiprocess_combined_list = multiprocess_combined_list
                               )
##    test.unpack_x(guess_x)
    test.optimize_x_IGC(guess = test.x[5])
    #test.get_mle(stringent_level = 'high')

    
if __name__ == '__main__':
##    parser = argparse.ArgumentParser()
##    parser.add_argument('--sim_num', required = True, help = 'simulation number')
##    parser.add_argument('--case', dest = 'case', required = True, help = 'which simulated data')
##    parser.add_argument('--guess', dest = 'guess', type=int, required = True, help = 'guess number')
##    
##    main(parser.parse_args())

    MyStruct = namedtuple('MyStruct', 'guess sim_num case')
    args = MyStruct(guess = 2, sim_num = 60, case = 'PSJSAnalyses')

    seq_file_list = np.loadtxt('../Harpak/missing_0_species_list.txt', dtype = str)[:1]
    pairs = [file_name.split(".")[:2] for file_name in seq_file_list]
    sim_num = args.sim_num
    case_folder = args.case

    alignment_file_list = ['./bootstrapReplicates/' + case_folder + '/' + '_'.join(paralog) \
                           + '/sim_' + str(sim_num) + '/' + '_'.join(paralog) + '_sim_' + str(sim_num) + '.fasta' for paralog in pairs]
    seq_index_file_list = ['../Harpak/prepared_input/' + seq_file.replace('.pos.seq.formatted', '').replace('.', '_') \
                           + '_seq_index.txt' for seq_file in seq_file_list]


    gene_to_orlg_file = '../Harpak/GeneToOrlg.txt'
    save_file = './save/' + case_folder + '/JS_Grand_save_sim_' + str(sim_num) + '.txt'

    tree_newick = '../Harpak/HarpakTree.newick'
    DupLosList = '../Harpak/HarpakDupLost.txt'
    terminal_node_list = ['Human', 'Chimp', 'Goril', 'Orang', 'Macaq', 'Mouse']
    node_to_pos = {'D1':0}

    ###### Now get HKY+PSJS-IGC estimates
    IGC_pm = 'One rate'
    pm_model = 'HKY'
    for paralog in pairs:
        log_paralog_folder = './log/' + case_folder + '/' + '_'.join(paralog)
        log_sim_folder = log_paralog_folder + '/sim_' + str(sim_num)
        save_paralog_folder = './save/' + case_folder + '/' + '_'.join(paralog)
        save_sim_folder = save_paralog_folder + '/sim_' + str(sim_num)
        summary_paralog_folder = './summary/' + case_folder + '/' + '_'.join(paralog)
        summary_sim_folder = summary_paralog_folder + '/sim_' + str(sim_num)

        if not os.path.isdir(log_paralog_folder):
            os.mkdir(log_paralog_folder)
        if not os.path.isdir(log_sim_folder):
            os.mkdir(log_sim_folder)
        if not os.path.isdir(save_paralog_folder):
            os.mkdir(save_paralog_folder)
        if not os.path.isdir(save_sim_folder):
            os.mkdir(save_sim_folder)
        if not os.path.isdir(summary_paralog_folder):
            os.mkdir(summary_paralog_folder)
        if not os.path.isdir(summary_sim_folder):
            os.mkdir(summary_sim_folder)

    save_file_list = ['./save/' + case_folder + '/' + '_'.join(paralog) + '/sim_' + str(sim_num) + '/JS_' + '_'.join(paralog) \
                           + '_guess_' + str(args.guess) + '_save.txt' for paralog in pairs]
    log_file_list = ['./log/' + case_folder + '/' + '_'.join(paralog) + '/sim_' + str(sim_num) + '/JS_' + '_'.join(paralog) \
                           + '_guess_' + str(args.guess) + '_log.txt' for paralog in pairs]
    summary_file_list = ['./summary/' + case_folder + '/' + '_'.join(paralog) + '/sim_' + str(sim_num) + '/JS_' + '_'.join(paralog) \
                           + '_guess_' + str(args.guess) + '_summary.txt' for paralog in pairs]

    rate_variation = False


    x_js = np.log([0.3,0.4,0.5,4.0, 3.0])
      
    test = Run_JS_Harpak_all(alignment_file_list, gene_to_orlg_file,
                 seq_index_file_list, 
                 tree_newick, DupLosList,
                 x_js, pm_model, IGC_pm, rate_variation,
                 node_to_pos, terminal_node_list,
                 save_file_list, log_file_list, save_file, 10)
##    test.get_mle(stringent_level = 'high')

    save_file = './save/' + case_folder + '/PSJS_Grand_save_sim_' + str(sim_num) + '.txt'

    save_file_list = ['./save/' + case_folder + '/' + '_'.join(paralog) + '/sim_' + str(sim_num) + '/PSJS_' + '_'.join(paralog) \
                           + '_guess_' + str(args.guess) + '_save.txt' for paralog in pairs]
    log_file_list = ['./log/' + case_folder + '/' + '_'.join(paralog) + '/sim_' + str(sim_num) + '/PSJS_' + '_'.join(paralog) \
                           + '_guess_' + str(args.guess) + '_log.txt' for paralog in pairs]
    summary_file_list = ['./summary/' + case_folder + '/' + '_'.join(paralog) + '/sim_' + str(sim_num) + '/PSJS_' + '_'.join(paralog) \
                           + '_guess_' + str(args.guess) + '_summary.txt' for paralog in pairs]
                           

    initial_tract_length_list = np.log([30.0, 200.0, 500.0])
    guess_lnp = -initial_tract_length_list[args.guess - 1]

    rate_variation = False
    true_x = np.loadtxt('./trueValues/PSJS_HKY_Harpak_One_rate_nonclock_save.txt')
    x_js = np.concatenate((true_x[:4], [true_x[4] - true_x[5] + guess_lnp, guess_lnp]))
    guess_x = np.concatenate((x_js, true_x[6:])) 


    # multiprocess_combined_list = [[0, 9, 12], [7, 5], [6, 1], [18, 8],
    #                               [2], [3], [4], [10], [11],
    #                               [13], [14], [15], [16], [17],
    #                               [19]]

    multiprocess_combined_list = None
      
    test = Run_PSJS_Harpak_all(alignment_file_list, gene_to_orlg_file,
                               seq_index_file_list,
                               tree_newick, DupLosList,
                               x_js, pm_model, IGC_pm, rate_variation,
                               node_to_pos, terminal_node_list,
                               save_file_list, log_file_list, save_file,
                               multiprocess_combined_list = multiprocess_combined_list
                               )
    if not os.path.isfile(test.save_file):
        test.unpack_x(guess_x)
    test.optimize_x_IGC(guess = test.x[5])
                         
