from HMMTract import *

if __name__ == '__main__':
    state_list = ['No IGC event (Si = 0)','At least one IGC event (Si > 0)']
    IGC_sitewise_lnL_file = './YLR406C_YDL075W_sitewise_lnL.txt'
    Force_sitewise_lnL_file = './Force_YLR406C_YDL075W_sitewise_lnL.txt'

    test = HMMTract(IGC_sitewise_lnL_file, Force_sitewise_lnL_file, state_list)
