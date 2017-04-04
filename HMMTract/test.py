from HMMTract import *

if __name__ == '__main__':
    state_list = ['No IGC event (Si = 0)','At least one IGC event (Si > 0)']
    IGC_sitewise_lnL_file = './YLR406C_YDL075W_sitewise_lnL.txt'
    Force_sitewise_lnL_file = './Force_YLR406C_YDL075W_sitewise_lnL.txt'
    Total_blen = 0.2

    test = HMMTract(IGC_sitewise_lnL_file, Force_sitewise_lnL_file, state_list,
                    Total_blen)

    self = test

    x = np.log([1.4, 0.4]) 
    lnL_array = test.Forward(True, x)
    result = test.get_mle()

    test.IGC_sitewise_lnL = np.array(test.IGC_sitewise_lnL) * 0.001

    lnL_array, state_array = test.Viterbi()
