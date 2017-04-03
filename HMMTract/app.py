## This code is copied from ST590C HW2
## It is modified to estimate IGC initiation rate and tract length by HMM
## Xiang Ji
## xji3@ncsu.edu
from numpy import *
import string
import copy
class HMM:
    def __init__ (self, Ptran, Emission, State_List, seq_dir, ds_dir):
        self.in_dir = seq_dir
        self.out_dir = ds_dir
        self.StateList = State_List # The ith state corresponds to the ith row
        self.Ptr = self.NormalizeRows(Ptran) # Pij is transition probability from state i to state j
        self.Emi = self.NormalizeRows(Emission)
        # columns: ACGT
        # ith row corresponds to the emission from ith state
        self.seq = self.ReadSeq()
        self.logL = 1.0
        self.threshold = 1e-200 # 1e-325 == 0.0 is True
        self.logthreshold = 500 # exp(-800)==0.0 is True
        self.logPtr = log(self.Ptr)
        self.logEmi = log(self.Emi)
        self.overflowshreshold = 1e300 # sys.floatinfo(max=1.7976931348623157e+308...

    def generateinitp(self, inits):
        assert(len(inits) == len(self.StateList))
        emit = self.seq[0]
        initp = ones(len(inits))
        emit_col = string.index('ACGY',emit)
        for i in range(len(self.StateList)):
            initp[i] = self.Emi[i,emit_col]*inits[i]

        return initp

    def ReadSeq(self):
        filename = self.in_dir + 'input.txt'
        with open(filename,'r') as f:
            seq = f.read()
        seq = string.replace(seq,'\n','')
        seq = string.replace(seq,' ','')
        seq = seq.upper()

        return seq

    def NormalizeRows(self,M):
        s = shape(M)
        for i in range(s[0]):
            M[i,:] = M[i,:]/sum(M[i,:])
            
        return M

        
    def Forward_step(self, emit, logp1):
        assert(len(logp1) == len(self.StateList))
        mi = logp1.max()
        #mi_idx = logp1.argmax()
        logp2 = ones(len(logp1))
        emit_col = string.index('ACGT',emit)
        
        for i in range(len(self.StateList)):
            logp2[i] = mi + self.logEmi[i,emit_col] + log(exp(-mi+logp1)*self.Ptr[:,i])

        return logp2

    def Forward_initial_jump(self, emit, p1):
        assert(len(p1)==len(self.StateList))
        p2 = ones(len(p1))
        emit_col = string.index('ACGT',emit)

        for i in range(len(self.StateList)):
            p2[i] = self.Emi[i,emit_col]*(p1*self.Ptr[:,i])

        return p2

    def Forward(self, inits):
        assert(len(inits) == len(self.StateList))
        initp = self.generateinitp(inits)
        logp = ones(len(self.StateList))
        if 0.0 in initp:
            if len(self.seq)>1:
                initp2 = self.Forward_initial_jump(self.seq[1],initp)
            else:
                return log(sum(initp))
            logp = log(initp2)
            for i in range(2,len(self.seq)):
                logp = self.Forward_step(self.seq[i],logp)
                if logp.max()-logp.min() >= self.logthreshold:
                    print 'There might be underflow error in the forward algorithm at step i = ', str(i)
                    print 'The logP array is ', logp

        else:
            logp = log(initp)
            for i in range(1,len(self.seq)):
                logp = self.Forward_step(self.seq[i],logp)
                if logp.max()-logp.min() >= self.logthreshold:
                    print 'There might be underflow error in the forward algorithm at step i = ', str(i)
                    print 'The logP array is ', logp
        mm = max(logp)
        result = mm + log(sum(exp(logp-mm)))

        return result

    def Viterbi_step(self, emit, logp1):
        assert(len(logp1) == len(self.StateList))
        logp2 = ones(len(logp1))
        path_state = []
        emit_col = string.index('ACGT',emit)

        for i in range(len(self.StateList)):
            allpossible = logp1 + self.logPtr[:,i].H + self.logEmi[i,emit_col]
            allpossible = squeeze(asarray(allpossible))
            logp2[i] = allpossible.max()
            path_state.append( [j for j,x in enumerate(allpossible) if x == allpossible.max()] )

        return logp2, path_state

    def Viterbi_initial_jump(self, emit, p1):
        assert(len(p1) == len(self.StateList))
        p2 = ones(len(p1))
        path_state = []
        emit_col = string.index('ACGT',emit)

        for i in range(len(self.StateList)):
            allpossible = multiply(p1 , self.Ptr[:,i].H) * self.Emi[i,emit_col]
            allpossible = squeeze(asarray(allpossible))
            p2[i] = allpossible.max()
            path_state.append( [j for j,x in enumerate(allpossible) if x == allpossible.max()] )

        return p2, path_state

    def Viterbi(self,inits):
        assert(len(inits) == len(self.StateList))
        initp = self.generateinitp(inits)
        logp = ones(len(self.StateList))
        path = []
        if 0.0 in initp:
            initp2, path_state = self.Viterbi_initial_jump(self.seq[1],initp)
            logp = log(initp2)
            path.append(path_state)
            for i in range(2, len(self.seq)):
                logp, path_state = self.Viterbi_step(self.seq[i],logp)
                path.append(path_state)
                if logp.min() <= -self.overflowshreshold: # since no exp() involved
                    print 'There might be underflow error in the viterbi algorithm at step i = ', str(i)

        else:
            logp = log(initp)
            for i in range(1, len(self.seq)):
                logp, path_state = self.Viterbi_step(self.seq[i],logp)
                path.append(path_state)
                if logp.min() <= -self.overflowshreshold: # since no exp() involved
                    print 'There might be underflow error in the viterbi algorithm at step i = ', str(i)

        maxlogp = logp.max()
        pathlastpos = [[i] for i,x in enumerate(logp) if x == logp.max()]
        allpath = self.tracebackpath(path,pathlastpos)

        return maxlogp, allpath

    def tracebackpath(self, path, pathlastpos):
        assert(len(path) == (len(self.seq)-1))
        allpath = copy.deepcopy(pathlastpos)
        for i in range(len(path)):
            allpath = self.traceback_step(len(self.seq)-1-i-1,path,allpath)

        return allpath


    def traceback_step(self, newpos, pathmatrix, inpath):
        assert((len(self.seq)-1)==len(pathmatrix))
        laststate = inpath[0][0]
        allonestate = True
        inpathnum = len(inpath)
        
        for i in range(1,inpathnum):
            allonestate = allonestate and inpath[i][0]
        assert(allonestate)

        newposstate = pathmatrix[newpos][laststate]
        outpathnum = inpathnum*len(newposstate)
        outpath = []
        for i in range(len(newposstate)):
           tmp = copy.deepcopy(inpath)
           [j.insert(0,newposstate[i]) for j in tmp]
           [outpath.append(j) for j in tmp]

        assert(len(outpath) == outpathnum)
        return outpath

    def printstatelist(self):
        for i in range(len(self.StateList)):
            print str(i), ' for state: ', self.StateList[i]

def displaypath(path):
    if len(path)==1:
        print 'There is only 1 optimal path'
    elif len(path)>1:
        print 'There are' + '\t' + str(len(path))+ '\t' + 'optimal paths'
        print
    else:
        print 'Please input a non-empty pathmatrix for func displaypath'
        print 
        
    for i in range(len(path)):
        print 'optimal path '+ '\t' + str(i+1) + '\t' + ':'
        print ''.join(str(e) for e in path[i])
        
        
        
        
        


        
        
                

    
        
        
        
        
                
                
        
            
        
