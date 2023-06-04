'''Quantum computation module
A module for quantum information threory 
Muhammad esmail hassani,
m s h j 6 7 (at) c h m a i l . i r
'''
from numpy import kron, array, matrix, ndarray, dot, random
from numpy.matrixlib.defmatrix import matrix as nmatrix
from numpy.random import choice
from sympy import simplify, Basic, log, exp, sqrt, I, E, pi
#from sympy import sec, sech, cos, cosh, sin, sinh
#from sympy import cot, coth, csc, csch, tan, tanh
#from sympy import asec, asech, acos, acosh, asin
#from sympy import asinh, acot, acoth, acsc, acsch, atan
from sympy import beta, floor, gamma, ln, atanh
import matplotlib
if matplotlib.get_backend()  ==  'agg':
    matplotlib.use('Qt4Agg')

from matplotlib import pyplot
from functools import reduce
import threading
import random
import re
e = E
zero = simplify(0)
one  = simplify(1)

class QDit:
    '''Help on class QDit in quantum module
'''
    __D1__ = {'1':[0,1],'0':[1,0],' - ':[0,1],' + ':[1,0]}

    
    __patternKetBell__ = '([|][01]+[>])'
    __patternKetX__    = '([|][+-]+[>])'

    __patternKetNat__ = '([|][0-9]+[>])'
    #__patternBraNat__ = '( < [0-9] + [|])'
    
    def __init__(self, state):
        '''Parameteres:
    Each of below structers:
        1. tuple('|0> + |3>  +  |7>  +  ...' , 5) #first  arg: state , secend arg: neglNumber
        2. tuple('2.1|2>  +  (1 + I)|5>  -  I|6>  +  ... + |101>' , 7)) #first  arg: state , secend arg: neglNumber
        3. 2 ** 2|00110101>  +  pi|00001100>
        4.  - I| +  +  -  +  -  +  +  + >
        5. log(2)|10111>  -  exp(3)|01011>
        5. sin(pi / 5)|1010>  +  atan(1.1)|1111>  +  exp(7.1)|0000>
        6. array([0,...,0,1,...,0,0,1,...,0,0]) # state of related matrix
        7. QDit!
'''
        
        if type(state) == str:
            state = state.replace(' ','')
            
            if re.search(self.__patternKetX__   ,state):
                self.__setstate__(state,self.__patternKetX__)
                self.dagger()
                self.__unicate__()
                self.state = Hadamard(self.len)(self).state
                

            elif re.search(self.__patternKetBell__,state):
                self.__setstate__(state,self.__patternKetBell__)
                self.dagger()

        elif isinstance(state,tuple) and len(state) == 2 and type(state[0]) == str and type(state[1]) == int:
            state = list(state)
            state[0] = state[0].replace(' ','')
            if re.search(self.__patternKetNat__,state[0]):
                self.__setstateNat__(state, self.__patternKetNat__)
                self.dagger()

            else:
                raise Exception("ket notation isn't in standard mod")

        elif isinstance(state,QDit):
            self.state = state.state

        elif isinstance(state,ndarray):
            try:
                self.state = array([complex(vi) for vi in state])
            except:
                self.state = array([complex(vi[0]) for vi in state])
            self.state = state
        else:
            raise Exception('QDit error')
        self.__unicate__()
        self.len = int(log(len(self.state),2))
        
    def __setstateNat__(self,state,patt):
        state = list(state)
        if state[0][0] == '|':
            state[0] = '1' + state[0]
        state[0] = state[0].replace('-|', '-1|')
        state[0] = state[0].replace('+|', '+1|')
        state,self.len = re.split(patt,state[0]),state[1]
        while '' in state: state.remove('')
        self.state = array([zero for i in range(2 ** self.len)])

        for i in range(1,len(state),2):
            if state[i - 1] in {'+','-'}:
                    state[i - 1] += '1'
            exec('v = ' + state[i - 1])
            self.state[int(state[i][1: - 1])] = simplify(locals()['v'])

    def __eq__(self,other):
        return (self.state == other.state).all()

    def __neq__(self,other):
        return (self.state != other.state).any()

    def __setstate__(self,state,patt):
        if state[0] == '|':
            state = '1' + state
        state = re.split(patt,state)
        while '' in state: state.remove('')
        self.len = len(state[1]) - 2
        
        self.state = array([zero for i in range(2 ** self.len)])
        for i in range(1,len(state),2):
            if state[i - 1] in {'+','-'}:
                state[i - 1] += '1'
            exec('v = ' + state[i - 1])
            vv = state[i][1: - 1].replace('+','0').replace('-','1')
            self.state[int(vv,2)] = simplify(locals()['v'])

    def __unicate__(self):
        d = zero
        for xi in self.state:
            d += abs(xi[0]) ** 2
        self.state /= sqrt(d)
        
    def dagger(self):
        '''return self dagger'''
        self.state = array([[vi.conjugate()] for vi in self.state])

    def __vdagger__(self):
        if self.KB == 'ket':
            return [vi[0].conjugate() for vi in self.state]
        return [[vi.conjugate()] for vi in self.state]

    def __Hadamard__(self):
        def t(V):
            return list(zip( * [list(Q.state)]))
        h = matrix([[1,1],[1, - 1]])
        l = int(log(len(self.state),2))
        H = matrix(reduce(kron,[h for i in range(l)]))
        
        self.state = (H * self.state).A

    def repr(self, rep = 'nat'):
        '''Parameters 
input:
    rep rep <-  {nat, bin, x}
Output:
    nat -> self representation by natural number
    bin -> self representation by {01}
     x   -> self representation by {+-}
'''
        if rep not in {'x', 'nat', 'bin'}:
            raise Exception("Representation mode must be one of {'x', 'nat', 'bin'}")
        
        represente = []
        if rep == 'x':
            preRep = [ vi[0] for vi in Hadamard(int(log(len(self.state),2)))(self)]
            j = 0
            while j < len(preRep):
                if preRep[j]:
                    x = bin(j)[2:]
                    x = ('0' * (int(log(len(preRep),2)) - len(x)) + x)
                    x = x.replace('0','+').replace('1','-')
                    represente += [str(preRep[j]) + '|' + x + '>']
                j += 1
            return ("QDit('" + '+'.join(represente) + "')").replace('+- ','-')

        ######################################

        represente = []
        if rep == 'nat':
            j = 0
            while j < len(self.state):
                if self.state[j][0]:
                    represente += [str(self.state[j][0]) + '|' + str(j) + '>']
                j += 1
            return ("QDit(('" + '+'.join(represente) + "' , "  + str(int(log(len(self.state), 2))) + "))").replace('+-','-')

        ##########################

        elif rep == 'bin':
            j = 0
            while j < len(self.state):
                if self.state[j][0]:
                    x = bin(j)[2:]
                    x = '0' * (int(log(len(self.state),2)) - len(x)) + x
                    represente += [str(self.state[j][0]) + '|' + x + '>']
                j += 1
            return ("QDit('" + '+'.join(represente) + "')").replace('+-','- ')

        '''elif 'part' in rep:

            if sum(args) != self.len:
                raise Exception('Sum of partitions length must be equal by number of quantum neglibles')

            L = [0]
            for ai in args:
                L += [L[ - 1] + ai]

            j = 0
            while j < len(self.state):
                if self.state[j][0]:
                    x = bin(j)[2:]
                    x = '0' * (self.len - len(x)) + x
                    x = ','.join([str(int(x[L[i - 1]:L[i]],2)) for i in range(1,len(L))])
                    represente += [str(self.state[j][0]) + '|' + x + '>']
                j += 1
            return ("QDit('" + '+'.join(represente) + "')").replace('+-','- ')'''

    def __repr__(self):
        return self.repr('nat')


    def __str__(self):
        return self.__repr__()



    def __mul__(self,other):
        if type(self) in {int,float,complex} or isinstance(self, Basic):
            return QDit(self * other.state)

        elif type(other) in {int,float,complex} or isinstance(other, Basic):
            return QDit(other * self.state)

        elif isinstance(self, QDit) and isinstance(other, QDit):
            return QDit(kron(self.state,other.state))

        elif isinstance(self, operator):
                return QDit(self.state * other)
        
        raise Exception('Bad multiplication requested')

    def __add__(self,other):
        return QDit(self.state + other.state)
    
    def __rmul__(self,other):
        return self.__mul__(other)

    def __radd__(self,other):
        return QDit(self.state + other.state)

    def __neg__(self):
        return self

    def __len__(self):
        return len(self.state)

    def __getitem__(self,index):
        return self.state[index]

    def dot(self,other):
        '''Parameters:
        QDit, QDit
    Return:
        float, internal product of input QDits stats vectors'''

        return dot(self.state.T[0], other.state.T[0])

class observer:
    '''Help on observer object in  quantum module'''
    def __init__(self, noiseBound  =  0):
        '''Parameters:
noiseBound
noise: an unwanted signal or a disturbance. in noise state vector, absolute state of each element is less than noiseBound.
if noiseBound  ==  0: observer dont imply any noise'''
        if noiseBound:
            self.__tip__  =  lambda x : noise(x,noiseBound)
        else:
            self.__tip__  =  lambda x : x

    def plot(self, inp, coordinate  =  'z'):
        '''plot:
Parameters:
    QDit, coordinate
    QDit: a QDit object
    coordinate: a subset of {'x', 'z', 'b'}
        'z'  <  -  standard coordinate that base  letters is {0, 1}
        'x'  <  -  oblique coordinate that base letters is { + ,  - }
        'b'  <  -  maximally entangled coordinate that base is {|i>  -  |(2  **  n  -  1) ^ i
Return:
   Plot chart of observ probeblities.'''
        coordinate  =  coordinate.lower()
        l  =  inp.len
        QDit  =  self.__tip__(inp)
        if (set(coordinate)  -  {'z','x','b'}):
            raise Exception('Coordinate have unexcepted charachter')
        if 'z' in coordinate:
            __ZProb__  = [abs(vi[0]) ** 2 for vi in inp.state]
            pyplot.plot(__ZProb__,label = 'Z base')

        if 'x' in coordinate:
            __XProb__ = [abs(vi[0]) ** 2 for vi in Hadamard(l)(inp).state]
            pyplot.plot(__XProb__,label = 'X base')

        if 'b' in coordinate:
            __BProb__ = [abs(vi[0]) ** 2 for vi in Bell(l)(inp).state]
            pyplot.plot(__BProb__,label = 'B base')

        pyplot.legend(loc = 'best')
        pyplot.show()        

    def output(self, inp, coordinate = 'Z'):
        '''output:
Parameters:
    QDit, coordinate
    QDit: a QDit object
    coordinate: a subset of {'x', 'z', 'b'}
        'z'  <  -  standard coordinate that base  letters is {0, 1}
        'x'  <  -  oblique coordinate that base letters is { + ,  - }
        'b'  <  -  maximally entangled coordinate that base is {a  -  b, a  +  b : for a, b in zBase (if a != b)}
Result:
    print observ detail in string mod
'''
        coordinate = coordinate.lower()
        l = inp.len
        inp = self.__tip__(inp)
        if (set(coordinate) - {'z','x','b'}):
            raise Exception('Coordinate have unexcepted charachter')
        ret = []
        if 'z' in coordinate:
            __ZProb__ = [abs(vi[0]) ** 2 for vi in inp.state]
            obs = choice(range(len(__ZProb__)),1,p = __ZProb__)[0]
            obs = obs,str(obs),bin(obs)[2:]
            obs  = obs[0],obs[1],'0' * (l - len(obs[2])) + obs[2]
            ret += ['']
            ret[ - 1] += "Observed state (Z base):\n\tQDit('|" + obs[2] + ">') = QDit(('|" + obs[1] + ">', " + str(l) + '))'
            ret[ - 1] += '\n\tPr(Obs = QDit(|' + obs[2] + '>))  =  ' + str(__ZProb__[obs[0]])
            
        if 'x' in coordinate:
            __XProb__ = [abs(vi[0]) ** 2 for vi in Hadamard(l)(inp).state]
            obs = choice(range(len(__XProb__)),1,p = __XProb__)[0]
            obs = obs,bin(obs)[2:].replace('0',' + ').replace('1',' - ')
            obs = obs[0],' + ' * (l - len(obs[1])) + obs[1]
            ret += ['']
            ret[ - 1] += "Observed state (X base):\n\tQDit('|" + obs[1] + ">')"
            ret[ - 1] += '\n\tPr(Obs = QDit(|' + obs[1] + '>))  =  ' + str(__XProb__[obs[0]])
            
        if 'b' in coordinate:
            BQ = Bell(l)(inp)
            __BProb__ = [abs(vi[0]) ** 2 for vi in BQ.state]
            obs = choice(range(len(__BProb__)),1,p = __BProb__)[0]
            ret += ['']
            ret[ - 1] += 'Observed state (B base):\n\tQDit(|B' + str(obs) + '>)'
            ret[ - 1] += '\n\tPr(Obs = |B' + str(obs) + '>)  =  ' + str(__BProb__[obs])

        print('\n\n'.join(ret))
        #return '\n\n'.join(ret)

class operator:
    '''Help on operator object in quantum module'''

    def __init__(self, value, n = 1):
        '''Parameters (valu, n):
       numpy.matrix, self tensor rank
   return
       An operator by state  **  n state'''
        if type(value) == operator:
            self.value = value.value
        else:
            if n == 1:
                self.value = matrix(value)
            elif n != 1:
                i = 0
                self.value = matrix([[1]])
                while i < n:
                    self.value = kron(self.value, value)
                    i += 1
        self.deg = int(log(len(self.value),2))

    def __len__(self):
        return self.deg

    def expand(self):
        '''Expand symbols that used in operator definition
O = operator([[e*(e+1) , 0],
                       [      0     , 1]])
O.expand = operator([[e + exp(2), 0],
                                   [      0        , 1]], dtype=object)
'''
        for i in range(len(self.value)):
            for j in range(len(self.value)):
                try:
                    self.value.A[i][j] = self.value.A[i][j].expand()
                except:
                    pass
    def __mul__(self,other):
        if isinstance(self, operator) and isinstance(other, operator):
            return operator(kron(self.value, other.value))
        
        elif isinstance(self, operator) and isinstance(other, QDit):
            return QDit(array(self.value * other.value))

        else:
            raise Exception('Bad multiplication requested')
    def __rmul__(self,other):
        return self.__mul__(other)

    def __pow__(self,n):
        i = 0
        ret = [[1]]
        while i < n:
            ret = kron(ret,self.value)
            i += 1
        return operator(ret)

    def __rpow__(self,other):
        return self.__pow__(other)

    def __getitem__(self,index):
        return self.value[index]

    def __add__(self,other):
        return operator(self.value + other.value)
    
    def __radd__(self,other):
        return self.__add__(other)

    def __repr__(self):
        return 'operator(' + repr(self.value)[7: - 1] + ')'

    def __str__(self):
        return 'operator(' + str(self.value) + ')'

    def __call__(self,inp):
        if isinstance(inp, QDit):
            return QDit(array(self.value * inp.state))
        if isinstance(inp, int):
            return self ** inp

def __sDiagonal__(A, B):
    ret, la, lb = [], len(A), len(B)
    for ai in A.A:
        ret += [list(ai) + [zero for i in range(lb)]]
    for bi in B.A:
        ret += [[zero for i in range(la)] + list(bi)]
    return ret

def __rDiagonal__(A,B):
    ret, la, lb = [], len(A[0]), len(B[0])
    for ai in A:
        ret += [[zero for i in range(lb)] + list(ai)]
    for bi in B:
        ret += [list(bi) + [zero for i in range(la)]]
    return ret

class Bell(operator):
    '''change to maximally Entangled bases, {a + b, a - b : for a , b in zBase (if a != b)}'''
    def __init__(self,degree):
        A = [[one * ( i == j) for i in range(2 ** degree)] for j in range(2 ** degree)]
        for i in range(2 ** (degree - 1), 2 ** degree):
            A[i][i] = -one
        A = operator(A)
        B = Not(degree)
        value = A + B
        super(Bell, self).__init__(value)

class Hadamard(operator):
    '''change Z base to X base'''
    def __init__(self,degree):
        super(Hadamard, self).__init__([[one, one],[one,  - one]],degree)

class QFT(operator):
    '''Quantum furiere transform'''
    def __init__(self,degree):
        d = 2 * pi * I / (2 ** degree)
        w = complex(e ** d)
        value = [[(w ** (i * j)) for i in range(2 ** degree)] for j in range(2 ** degree)]
        super(QFT, self).__init__(value)

class Identity(operator):
    '''Identity operator'''
    def __init__(self,degree):
        super(Identity, self).__init__([[one, zero],[zero, one]],degree)
    
class   Not(operator):
    '''Flipp all neglible'''
    def __init__(self,degree):
        super(Not, self).__init__([[zero, one],[one, zero]],degree)

class sqrNot(operator):
    def __init__(self):
        super(sqrNot, self).__init__([[one + I,one - I],[one - I,one + I]])

class CNot(operator):
    def __init__(self):
        value = __sDiagonal__(Identity(1).state, Not(1).state)
        super(CNot, self).__init__(value)

class pauliy(operator):
    def __init__(self,degree):
        super(Pauliy, self).__init__([[zero,  - I],[I, zero]],degree)

class pauliz(operator):
    def __init__(self,degree):
        super(Pauliz, self).__init__([[one, zero],[zero,  - one]],degree)


class phaseShift(operator):
    def __init__(self,degree,phi):
        r = e ** (I * phi)
        super(PhaseShift, self).__init__([[one, zero],[zero, r]],degree)

class swap(operator):
    def __init__(self):
        value = [[one, zero, zero, zero], [zero, zero, one, zero],
               [zero, one, zero, zero], [zero, zero, zero, one]]
        super(swap, self).__init__(value)

class squarSwap(operator):
    def __init__(self):
        value = [[one, zero, zero, zero],[zero, 0.5 * (1 + I), 0.5 * (1 - I), zero],
               [zero, 0.5 * (1 - I), 0.5 * (1 + I), zero],[zero, zero, zero, one]]
        super(SquarSwap, self).__init__(value)

class toffali(operator):
    def __init__(self,degree):
        if degree < 3:
            raise Exception('Toffali.degree > =  3')
        value = __sDiagonal__(Identity(degree - 1).value, Not(1).value)
        super(toffali, self).__init__(value)


class fredkin(operator):
    def __init__(self,degree):
        value = Identity(degree).value
        value[-2, -2] = value[-3, -3]= 0
        value[-2, -3] = value[-3, -2] = 1
        super(fredkin, self).__init__(value)


class ising(operator):
    def __init__(self,phi):
        value = [[one, zero, zero,  - I * e ** (I * phi)],[zero, one, - I,zero],
               [zero , - I, one, zero],[ - I * e ** (I * phi), zero, zero, one]]
        super(Ising, self).__init__(value)


def noise(inp, noiseBound = 0.1):
    '''Parematers:
    QDite, norseBound
Return:
    QDite'''
    degree = len(inp)
    rand = random.rand(degree, degree) * noiseBound
    return operator(rand)(inp)

############################################
############################################
import sys

def get_size(obj, seen=None):
    """Recursively finds size of objects"""
    size = sys.getsizeof(obj)
    if seen is None:
        seen = set()
    obj_id = id(obj)
    if obj_id in seen:
        return 0
    # Important mark as seen *before* entering recursion to gracefully handle
    # self-referential objects
    seen.add(obj_id)

    try:
        size += get_size(obj.state, seen)
    except Exception as a:
        pass

    if isinstance(obj, dict):
        size += sum([get_size(v, seen) for v in obj.states()])
        size += sum([get_size(k, seen) for k in obj.keys()])
    elif hasattr(obj, '__dict__'):
        size += get_size(obj.__dict__, seen)
    elif hasattr(obj, '__iter__') and not isinstance(obj, (str, bytes, bytearray)):
        size += sum([get_size(i, seen) for i in obj])
    #elif isinstance(obj, ndarray):
    try:
        for oi in obj:
            size += get_size(oi, seen)
    except:
        pass
    return size

