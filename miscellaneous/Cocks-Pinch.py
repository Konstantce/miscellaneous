import math

EMBD_DEGREES = [6, 12, 24]
START_DISCR = -3
STATIC_R = 21888242871839275222246405745257275088548364400416034343698204186575808495617
#STATIC_R = 982451653
STATIC_R = 13381

def check_embedding_degree(r):
    for k in EMBD_DEGREES:
        if (r - 1) % k == 0:
            return True
    return False


def get_root_of_unity(r, k):
    field = GF(r)
    gen = field.multiplicative_generator()
    a = (r - 1) / (k)
    root_of_unity = gen ^ a
    return root_of_unity


def check_if_fundamental_Disc(d):
    if -d % 16 not in [3, 4, 7, 8, 11, 15]:
        return False
    M = squarefree_part(d)
    L = d / M
    while L % 2 == 0:
        L = L / 2
    return True if L == 1 else False
    
def get_params(r):
    D = START_DISCR
    field = GF(r)
    
    while True:
        if not check_if_fundamental_Disc(D):
            D -= 1
            continue
            
        if kronecker(D, r) != 1:
            D -= 1
            continue
        
        for k in EMBD_DEGREES:
            if (r - 1) % k == 0:
                g = get_root_of_unity(r, k)
                t_ = g + 1
                u_ = (t_ - 2) / field(D).sqrt()
    
                t = int(t_)
                u = int(u_)
            
                p = (t^2 - D * u^2) / 4

                if p in ZZ and p in Primes():
                    return D, p, int(g), k
        D -=1

        
def find_nonresidue(field):
    p = field.characteristic()
   
    while True:
        x = field.random_element()
        if x.is_square():
            continue
        if ( p % 3 == 1) and (x^((p-1) / 3) == field(1)):
            continue
        
        return x


def get_curve_params(D, p):
    #return all possible pairs!
    
    field = GF(p)
    g = find_nonresidue(field)
    
    if D == -4:
        res = []
        for k in xrange(4):
            pair = (int(-g^k), 0)
            res.append(pair)
        return res

    if D == -3:
        res = []
        for k in xrange(6):
            pair = (0, int(-g^k))
            res.append(pair)
        return res
    

    Hilbert_poly = hilbert_class_polynomial(D)
    #reduced_poly = Hilbert_poly.change_ring(field)
    #j = reduced_poly.roots()[0][0]
    j = Hilbert_poly.any_root(field)
    
    c = j/(j - field(1728))
    r = field(-3)*c
    s = field(2)*c
    
    return [(int(r), int(s)), (int(r * g^2), int(s * g^3))]


def Cocks_Pinch(r):
    if not check_embedding_degree(r):
        print "Unsatisfiable"
        return False

    D, p, g, k = get_params(r)
    coeffs_list = get_curve_params(D, p)
    
    return p, coeffs_list, k, D


def test_ring_change():
    c = PolynomialRing(ZZ, 'w')
    poly = c("w - 5")
    reduced_poly = poly.change_ring(GF(3))
    print reduced_poly.roots()
    

def select_curve(p, coeffs_list, k, r):
    field = GF(p)
    for (A, B) in coeffs_list:
        curve = EllipticCurve(field, [A, B])
        order = curve.cardinality()
        if order % r == 0:
            
            print "Curve found!"
            return A, B
        
    
def find_group_generators(p, k, A, B, r):
    K.<a> = GF(p^k)
    E = EllipticCurve(GF(p), [A, B])
    EK = E.base_extend(K)
    
    cofactor = int(EK.cardinality() / r)
    
    G1 = E.random_point()
    while G1.order() % r != 0:
        G1 = E.random_point()
    G1 = G1 * int(G1.order() / r)

    if G1.order() != r:
        print "G1 Error!"
        return False
    
    G2 = EK.random_point()
    while G2.order() % r != 0:
        G2 = EK.random_point()
    G2 = G2 * int(G2.order() / r)

    if G2.order() != r:
        print "G1 Error!"
        return False
    
    P = EK(G1)
    Gt = P.tate_pairing(G2, r, k)
    
    if Gt.multiplicative_order() != r:
        print "Gt Error!"
        return False

    return G1, G2, Gt
    
    
def print_params(p, A, B, k, G1, G2, Gt):
    print "Curve base field: ", p
    print "A: ", A
    print "B: ", B
    print "embedding degree: ", k
    print "generator of first source group: ", G1
    print "generator of second source group: ", G2
    print "generator of destination group: ", Gt
    
        
    
if __name__ == "__main__":   
    
    p, coeffs_list, k, D = Cocks_Pinch(STATIC_R)
    A, B = select_curve(p, coeffs_list, k, STATIC_R)
    G1, G2, Gt = find_group_generators(p, k, A, B, STATIC_R)
    
    print_params(p, A, B, k, G1, G2, Gt)
