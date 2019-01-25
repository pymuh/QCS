'''Quantum computation module
A module for quantum information threory 
Muhammad esmail hassani,
m s h j 6 7 (at) c h m a i l . i r
'''
from numpy import kron, array, matrix, ndarray, dot
from numpy.matrixlib.defmatrix import matrix as nmatrix
from numpy.random import choice
from sympy import simplify, log, exp, sqrt, I, E, pi
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
    
    def __init__(self,value):
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
        if type(value) == str:
            value = value.replace(' ','')

            if re.search(self.__patternKetX__   ,value):
                self.__setValue__(value,self.__patternKetX__)
                self.dager()
                self.__unicate__()
                self.value = Hadamard(self.len)(self).value
                

            elif re.search(self.__patternKetBell__,value):
                self.__setValue__(value,self.__patternKetBell__)
                self.dager()

        elif isinstance(value,tuple) and len(value) == 2 and type(value[0]) == str and type(value[1]) == int:
            
            if re.search(self.__patternKetNat__,value[0]):
                self.__setValueNat__(value,self.__patternKetNat__)
                self.dager()

            else:
                raise Exception("ket notation isn't in standard mod")

        elif isinstance(value,QDit):
            self.value = value.value

        elif isinstance(value,ndarray):
            try:
                self.value = array([complex(vi) for vi in value])
            except:
                self.value = array([complex(vi[0]) for vi in value])
            self.value = value
        else:
            raise Exception('QDit error')
        self.__unicate__()
        self.len = int(log(len(self.value),2))
        
    def __setValueNat__(self,value,patt):
        value = list(value)
        if value[0][0] == '|':
            value[0] = '1' + value[0]
        value,self.len = re.split(patt,value[0]),value[1]

        while '' in value: value.remove('')
        self.value = array([zero for i in range(2 ** self.len)])

        for i in range(1,len(value),2):
            if value[i - 1] in {'+','-'}:
                    value[i - 1] += '1'
            exec('v = ' + value[i - 1])
            self.value[int(value[i][1: - 1])] = simplify(locals()['v'])

    def __eq__(self,other):
        return (self.value == other.value).all()

    def __neq__(self,other):
        return (self.value != other.value).any()

    def __setValue__(self,value,patt):
        if value[0] == '|':
            value = '1' + value
        value = re.split(patt,value)
        while '' in value: value.remove('')
        self.len = len(value[1]) - 2
        
        self.value = array([zero for i in range(2 ** self.len)])
        for i in range(1,len(value),2):
            if value[i - 1] in {'+','-'}:
                value[i - 1] += '1'
            exec('v = ' + value[i - 1])
            vv = value[i][1: - 1].replace('+','0').replace('-','1')
            self.value[int(vv,2)] = simplify(locals()['v'])

    def __unicate__(self):
        d = zero
        for xi in self.value:
            d += abs(xi[0]) ** 2
        self.value /= sqrt(d)
        
    def dager(self):
        '''return self dager'''
        self.value = array([[vi.conjugate()] for vi in self.value])

    def __vDager__(self):
        if self.KB == 'ket':
            return [vi[0].conjugate() for vi in self.value]
        return [[vi.conjugate()] for vi in self.value]

    def __Hadamard__(self):
        def t(V):
            return list(zip( * [list(Q.value)]))
        h = matrix([[1,1],[1, - 1]])
        l = int(log(len(self.value),2))
        H = matrix(reduce(kron,[h for i in range(l)]))
        
        self.value = (H * self.value).A

    def repr(self,rep = 'nat'):
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
            preRep = [ vi[0] for vi in Hadamard(int(log(len(self.value),2)))(self)]
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
            while j < len(self.value):
                if self.value[j][0]:
                    represente += [str(self.value[j][0]) + '|' + str(j) + '>']
                j += 1
            return ("QDit(('" + '+'.join(represente) + "' , "  + str(int(log(len(self.value), 2))) + "))").replace('+-','-')

        ##########################

        elif rep == 'bin':
            j = 0
            while j < len(self.value):
                if self.value[j][0]:
                    x = bin(j)[2:]
                    x = '0' * (int(log(len(self.value),2)) - len(x)) + x
                    represente += [str(self.value[j][0]) + '|' + x + '>']
                j += 1
            return ("QDit('" + '+'.join(represente) + "')").replace('+-','- ')

        '''elif 'part' in rep:

            if sum(args) != self.len:
                raise Exception('Sum of partitions length must be equal by number of quantum neglibles')

            L = [0]
            for ai in args:
                L += [L[ - 1] + ai]

            j = 0
            while j < len(self.value):
                if self.value[j][0]:
                    x = bin(j)[2:]
                    x = '0' * (self.len - len(x)) + x
                    x = ','.join([str(int(x[L[i - 1]:L[i]],2)) for i in range(1,len(L))])
                    represente += [str(self.value[j][0]) + '|' + x + '>']
                j += 1
            return ("QDit('" + '+'.join(represente) + "')").replace('+-','- ')'''

    def __repr__(self):
        return self.repr('nat')


    def __str__(self):
        return self.__repr__()


    def str(self,C = 'B'):
        '''Return: str, self string represantation'''
        if C == 'B':
            return str(self)
        if C == 'X':
            return str(Hadamard(self))


    def __mul__(self,other):
        
        if type(self) in {int,float,complex}:
            return QDit(self * other.value)

        elif type(other) in {int,float,complex}:
            return QDit(other * self.value)

        elif isinstance(self, QDit) and isinstance(other, QDit):
            return QDit(kron(self.value,other.value))

        elif type(self) == operator:
                return QDit(self.value * other)
        
        raise Exception('Bad multiplication requested')

    def __add__(self,other):
        return QDit(self.value + other.value)
    
    #def __rmul__(self,other):
        #return self.value * other.value

    def __radd__(self,other):
        return QDit(self.value + other.value)

    def __neg__(self):
        return self

    def __len__(self):
        return len(self.value)

    def __getitem__(self,index):
        return self.value[index]

    def dot(self,other):
        '''Parameters:
        QDit, QDit
    Return:
        float, internal product of input QDits stats vectors'''

        if self.KB == 'bra' and other.B == 'ket':
            return self * other

class observer:
    '''Help on observer object in  quantum module'''
    def __init__(self, noiseBound  =  0):
        '''Parameters:
noiseBound
noise: an unwanted signal or a disturbance. in noise state vector, absolute value of each element is less than noiseBound.
if noiseBound  ==  0: observer dont imply any noise'''
        if noiseBound:
            self.__tip__  =  lambda x : self.__noise__(x,noiseBound)
        else:
            self.__tip__  =  lambda x : x

    def __noise__(self, QDit, noiseBound = 0.1):
        degree  =  len(QDit)
        rand  =  lambda :simplify(random.uniform( - noiseBound,  + noiseBound))
        return QDit(array([[simplify(float(vi[0]  +  rand()))] for vi in QDit]))

    def plot(self, QDit, coordinate  =  'z'):
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
        l  =  int(log(len(QDit.value),2))
        QDit  =  self.__tip__(QDit)
        if (set(coordinate)  -  {'z','x','b'}):
            raise Exception('Coordinate have unexcepted charachter')
        if 'z' in coordinate:
            __ZProb__  = [abs(vi[0]) ** 2 for vi in QDit.value]
            pyplot.plot(__ZProb__,label = 'Z base')

        if 'x' in coordinate:
            __XProb__ = [abs(vi[0]) ** 2 for vi in Hadamard(l)(QDit).value]
            pyplot.plot(__XProb__,label = 'X base')

        if 'b' in coordinate:
            __BProb__ = [abs(vi[0]) ** 2 for vi in Bell(l)(QDit).value]
            pyplot.plot(__BProb__,label = 'B base')

        pyplot.legend(loc = 'best')
        pyplot.show()        

    def output(self,QDit, coordinate = 'Z'):
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
        l = int(log(len(QDit.value),2))
        QDit = self.__tip__(QDit)
        if (set(coordinate) - {'z','x','b'}):
            raise Exception('Coordinate have unexcepted charachter')
        ret = []
        if 'z' in coordinate:
            __ZProb__ = [abs(vi[0]) ** 2 for vi in QDit.value]
            obs = choice(range(len(__ZProb__)),1,p = __ZProb__)[0]
            obs = obs,str(obs),bin(obs)[2:]
            obs  = obs[0],obs[1],'0' * (l - len(obs[2])) + obs[2]
            ret += ['']
            ret[ - 1] += "Observed state (Z base):\n\tQDit('|" + obs[2] + ">') = QDit(('|" + obs[1] + ">', " + str(l) + '))'
            ret[ - 1] += '\n\tPr(Obs = QDit(|' + obs[2] + '>))  =  ' + str(__ZProb__[obs[0]])
            
        if 'x' in coordinate:
            __XProb__ = [abs(vi[0]) ** 2 for vi in Hadamard(l)(QDit).value]
            obs = choice(range(len(__XProb__)),1,p = __XProb__)[0]
            obs = obs,bin(obs)[2:].replace('0',' + ').replace('1',' - ')
            obs = obs[0],' + ' * (l - len(obs[1])) + obs[1]
            ret += ['']
            ret[ - 1] += "Observed state (X base):\n\tQDit('|" + obs[1] + ">')"
            ret[ - 1] += '\n\tPr(Obs = QDit(|' + obs[1] + '>))  =  ' + str(__XProb__[obs[0]])
            
        if 'b' in coordinate:
            BQ = Bell(l)(QDit)
            __BProb__ = [abs(vi[0]) ** 2 for vi in BQ.value]
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
       An operator by value  **  n state'''
        if type(value) == operator:
            self.value = value.value
        else:
            if n == 1:
                self.value = matrix(value)
            elif n != 1:
                i = 0
                self.value = matrix([[1]])
                while i < n:
                    self.value = kron(self.value,value)
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
            return operator(kron(self.value,other.value))
        
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
            return QDit(array(self.value * inp.value))
        if isinstance(inp, int):
            return self ** inp

def __sDiagonal__(A,B):
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
        super(Bell,self).__init__(value)

class Hadamard(operator):
    '''change Z base to X base'''
    def __init__(self,degree):
        super(Hadamard,self).__init__([[one, one],[one,  - one]],degree)

class QFT(operator):
    '''Quantum furiere transform'''
    def __init__(self,degree):
        d = 2 * pi * I / (2 ** degree)
        w = complex(e ** d)
        value = [[(w ** (i * j)) for i in range(2 ** degree)] for j in range(2 ** degree)]
        super(QFT,self).__init__(value)

class Identity(operator):
    '''Identity operator'''
    def __init__(self,degree):
        super(Identity,self).__init__([[one, zero],[zero, one]],degree)
    
class Not(operator):
    '''Flipp all neglible'''
    def __init__(self,degree):
        super(Not,self).__init__([[zero, one],[one, zero]],degree)

class SqrNot(operator):
    def __init__(self):
        super(SqrNot,self).__init__([[one + I,one - I],[one - I,one + I]])

class CNot(operator):
    def __init__(self):
        value = __sDiagonal__(Identity(1).value.A,Not(1).value.A)
        super(CNot,self).__init__(value)

class Pauliy(operator):
    def __init__(self,degree):
        super(Pauliy,self).__init__([[zero,  - I],[I, zero]],degree)

class Pauliz(operator):
    def __init__(self,degree):
        super(Pauliz,self).__init__([[one, zero],[zero,  - one]],degree)


class PhaseShift(operator):
    def __init__(self,degree,phi):
        r = e ** (I * phi)
        super(PhaseShift,self).__init__([[one, zero],[zero, r]],degree)

class Swap(operator):
    def __init__(self):
        value = [[one, zero, zero, zero], [zero, zero, one, zero],
               [zero, one, zero, zero], [zero, zero, zero, one]]
        super(Swap,self).__init__(value)

class SquarSwap(operator):
    def __init__(self):
        value = [[one, zero, zero, zero],[zero, 0.5 * (1 + I), 0.5 * (1 - I), zero],
               [zero, 0.5 * (1 - I), 0.5 * (1 + 1j), zero],[zero, zero, zero, one]]
        super(SquarSwap,self).__init__(value)

class Toffali(operator):
    def __init__(self,degree):
        if degree < 3:
            raise Exception('Toffali.degree > =  3')
        value = __sDiagonal__(Identity(degree - 1).value,Not(1).value)
        super(Toffali,self).__init__(value)


class Fredkin(operator):
    def __init__(self,degree):        
        value = __sDiagonal__(Identity(degree - 2).value,Swap().value)
        super(Fredkin,self).__init__(value)


class Ising(operator):
    def __init__(self,phi):
        value = [[one, zero, zero,  - I * e ** (I * phi)],[zero, one, - I,zero],
               [zero , - I, one, zero],[ - I * e ** (I * phi), zero, zero, one]]
        super(Ising,self).__init__(value)


def noise(QDit,noiseBound = 0.1):
    '''Parematers:
    QDite, norseBound
Return:
    QDite'''
    degree = len(QDit)
    rand = lambda :simplify(random.uniform( - noiseBound, + noiseBound))
    return QDit(array([[simplify(float(vi[0] + rand()))] for vi in QDit]))

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
        size += get_size(obj.value, seen)
    except Exception as a:
        pass

    if isinstance(obj, dict):
        size += sum([get_size(v, seen) for v in obj.values()])
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

T = []
from time import time
#a = QDit(('1|0>' , 1))
#print(get_size(a))
a = QDit('|' + '0' * 4 + '>')
'''for i in range(1,30):
        t = time()
        a = QDit('|' + '0' * i + '>')
        t = time() - t
        T += [[i, round(t, 5)]]
        print(T[-1])
        del a
'''
