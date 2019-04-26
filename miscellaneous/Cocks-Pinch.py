import math

EMBD_DEGREES = [6, 12, 24]
START_DISCR = -3
#STATIC_R = 21888242871839275222246405745257275088548364400416034343698204186575808495617
STATIC_R = 982451653

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

    
def get_params(r):
    D = START_DISCR
    field = GF(r)
    
    while True:
        if D % 4 == 2 or D % 4 == 3:
            D -= 1
            continue
            
        if kronecker(-D, r) != 1:
            D -= 1
            continue
        
        for k in EMBD_DEGREES:
            if (r - 1) % k == 0:
                g = get_root_of_unity(r, k)
                t_ = g + 1
                u_ = (t_ - 2) / field(-D).sqrt()
    
                t = int(t_)
                u = int(u_)
                p = (t^2 + D * u^2) / 4

                if p in ZZ and p in Primes():
                    print D
                    return D, p, g, k
        D -=1
        print D

    
def get_curve_params(D, p):
    if D == -4:
        return (-1, 0)
    if D== -3:
        return (0, -1)
    
    field = GF(p)
    Hilbert_poly = hilbert_class_polynomial(D)
    reduced_poly = Hilbert_poly.change_ring(field)
    j = reduced_poly.roots()[0][0]
    print j
    
    c = j/(j - field(1728))
    r = int(field(-3)*c)
    s = int(field(2)*c)
    
    return (r, s)


def Cocks_Pinch(r):
    if not check_embedding_degree(r):
        print "Unsatisfiable"
        return False

    D, p, g, k = get_params(r)

    A, B = get_curve_params(D, p)
    
    base_field = GF(p)
    extension_field = GF(p^k, name = 't')
    ext_field_modulus = extension_field.modulus()
    
    return p, A, B, k


def test_ring_change():
    c = PolynomialRing(ZZ, 'w')
    poly = c("w - 5")
    reduced_poly = poly.change_ring(GF(3))
    print reduced_poly.roots()
    
    
def check_curve(p, A, B, k, l):
    field = GF(p)
    curve = EllipticCurve(field, A, B)
    print "HERE!"
    order = curve.cardinality()
    print order, (order % l)

    
if __name__ == "__main__":
    
    p, A, B, k = Cocks_Pinch(STATIC_R)
    check_curve(p, A, B, k, STATIC_R)
    
    print "Base Field = ", p
    print "A = ", A
    print "B = ", B
    print "k = ", k
