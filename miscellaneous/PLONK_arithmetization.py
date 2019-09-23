
               import itertools
import collections

#PLONK arithmetizer
p = 257
field = GF(p)
omega = field.multiplicative_generator()
poly_ring = PolynomialRing(field, 'x')

def next_power_of_2(x):  
    return 1 if x == 0 else 2**(x - 1).nbits()


def construct_interpolation_poly(domain, values):
    return poly_ring.lagrange_polynomial(itertools.izip_longest(domain, values, fillvalue=0))


class Permutation:
    def __init__(self, cycles):
        self.n = sum([len(x) for x in cycles], 0)
        self.func = {}
        for cycle in cycles:
            a = collections.deque(cycle)
            b = collections.deque(cycle)
            b.rotate(1)
            for (i, j) in zip(a, b):
                self.func[i] = j
    
    def __call__(self, i):
        if i < 0 or i >= self.n:
            raise Exception("Permutation argument is incorrect")
        return self.func[i]
    
    def size(self):
        return self.n


class Wire:
    id_count =  0
    
    def __init__(self, is_public = False):
        self.id = Wire.id_count
        self.is_public = is_public
        Wire.id_count += 1
        
    @staticmethod
    def get_num_of_wires():
        return id_count
    
    def get_id(self):
        return self.id
        
        
class Constraint:
    #general constraint
    def __init__(self, q_l, q_r, q_o, q_m, q_c, a_wire, b_wire, c_wire):
        self.q_l = field(q_l)
        self.q_r = field(q_r) 
        self.q_o = field(q_o)
        self.q_m = field(q_m)
        self.q_c = field(q_c)
        self.a_wire = a_wire
        self.b_wire = b_wire
        self.c_wire = c_wire
        

class ConstraintSystem:
    
    #l here means the total number of public inputs used by the constraint system
    def __init__(self, l):
        self.wires = []
        self.public_wires = []
        self.constraints = []
        self.total_num_public_wires = l
        self.cur_num_public_wires = 0
    
                         
    def add_wire(self, is_public = False):
        wire = Wire(is_public)
        self.wires.append(wire)
        if is_public:
            self.public_wires.append(wire)
            self.cur_num_public_wires += 1
                         
        return wire
    
                         
    def add_constraint(self, q_l, q_r, q_o, q_m, q_c, a_wire, b_wire, c_wire):
        if any(elem not in self.wires  for elem in (a_wire, b_wire, c_wire)):
            raise Exception("undefined wire is used in the constraint")
        
        cnstr = Constraint(q_l, q_r, q_o, q_m, q_c, a_wire, b_wire, c_wire)
        self.constraints.append(cnstr)
                         
        return cnstr
    
                         
    def check_if_prepared(self):
        if self.cur_num_public_wires != self.total_num_public_wires:
            raise Exception("Inconsistent number of wires")
                         
        if len(self.constraints) < self.total_num_public_wires:
            raise Exception("Not enough constraints!")
        
        temp_storage = []
                         
        for i in xrange(self.total_num_public_wires):
            c = self.constraints[i]
            if any([c.q_l != field(1), c.q_r != field(0), c.q_0 != field(0), c.q_m != field(0), c.q_c != field(0)]):
                   raise Exception("Selectors are incorrect")
                   
            
            if c.a_wire not in self.public_wires or c.a_wire in temp_storage:
                raise Exception("Wires are incorrect")
            temp_storage.append(c.a_wire)
            
                        
    def construct_poly_encoding(self):
        #check if constraint system is in prepared form
        self.check_if_prepared()
                   
        m = len(self.wires)
        n = len(self.constraints)
                   
        #we do want n+1 to be a power of 2
        #if it is not the case we increase n to the closest higher power of 2
        # additional constaints are of the same form: all of the selectors is 0 and wire is used is always the first
        
        n = next_power_of_2(n+1) - 1
                   
        if (p-1) % (n+1) != 0:
            raise Exception("used field is too small to embed all of the constraints")
        d = int((p-1)/ (n+1))           
        
        g = omega ** d
                   
        #if this asserion is really important:
        if (m < n) or (m > 2*n):
            raise Exception("unbalanced contraint system")
                   
        domain = [g**k for k in xrange(1, n+1)]
                   
        #construct selector polynomials: q_L, q_R, q_C, q_M, q_O
        q_L = construct_interpolation_poly(domain, [c.q_l for c in self.constraints])
        q_R = construct_interpolation_poly(domain, [c.q_r for c in self.constraints])
        q_M = construct_interpolation_poly(domain, [c.q_m for c in self.constraints])
        q_C = construct_interpolation_poly(domain, [c.q_c for c in self.constraints])
        q_O = construct_interpolation_poly(domain, [c.q_o for c in self.constraints])
                   
        #construct S_{id1}, S_{id2}, S_{id3}
        # S_{id1}(g^i) = i
        # S_{idj}(x) = S_{id1}(x) + (j − 1) · n      
        S1 = construct_interpolation_poly(domain, [i for i in xrange(1, n+1)])
        S2 = S1 + n
        S3 = S1 + 2*n
                   
        #construct permutation polynomials P1, P2, P3
        #we start with defining the permutation \delta    
        perm = {k: [] for k in xrange(m)}
        for i, cnstr in enumerate(self.constraints):
            perm[cnstr.a_wire.get_id()].append(i)
            perm[cnstr.b_wire.get_id()].append(i+n)
            perm[cnstr.c_wire.get_id()].append(i+2*n)
                   
        #add elements from dumb padding constraints
        nc = len(self.constraints)
        temp = n - nc
        if temp > 0:
            perm[0] += [i for i in xrange(nc, n)]
            perm[0] += [i for i in xrange(nc+n, 2*n)]
            perm[0] += [i for i in xrange(nc+2*n, 3*n)]
                   
        #construct permutation as a function
        perm = Permutation(list(perm.values()))
        if perm.size() != 3*n:
             raise Exception("Permutation size is incorrect")
        
        #P_{j}(g^i) = perm((j−1)·n+i)
        P1 = construct_interpolation_poly(domain, [perm(i) for i in xrange(n)])
        P2 = construct_interpolation_poly(domain, [perm(i+n) for i in xrange(n)])
        P3 = construct_interpolation_poly(domain, [perm(i+2*n) for i in xrange(n)])
                   
        return q_L, q_R, q_M, q_C, q_O, S1, S2, S3, P1, P2, P3
    
    
def print_polys(q_L, q_R, q_M, q_C, q_O, S1, S2, S3, P1, P2, P3):
    print "q_L: ", q_L
    print "q_R: ", q_R
    print "q_M: ", q_M
    print "q_C: ", q_C
    print "q_O: ", q_O
    print "S1: ", S1
    print "S2: ", S2
    print "S3: ", S3
    print "P1: ", P1
    print "P2: ", P2
    print "P3: ", P3
                               
        
def build_test_constraint_system():
    #all variables are private: a, b, c
    #2*a-b = 0
    #a*b = c
    
    system = ConstraintSystem(0)
    a = system.add_wire()
    b = system.add_wire()
    c = system.add_wire()
    
    system.add_constraint(2, -1, 0, 0, 0, a, b, c)
    system.add_constraint(0, 0, 1, 1, 0, a, b, c)
    
    polys = system.construct_poly_encoding()
    print_polys(*polys)
    
    
build_test_constraint_system()




               
