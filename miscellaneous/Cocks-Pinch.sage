import time
import sys
import json

    
def is_valid_curve(q,t,r,k,D): 
    """
    Description:
    
        Tests that (q,t,r,k,D) is a valid elliptic curve
    
    Input:
    
        q - size of prime field
        t - trace of Frobenius
        r - size of prime order subgroup
        k - embedding degree
        D - (negative) fundamental discriminant
    
    Output:
    
        bool - true iff there exists an elliptic curve over F_q with trace t, a subgroup of order r with embedding degree k, and fundamental discriminant D
    
    """
    if q == 0 or t == 0 or r == 0 or k == 0 or D == 0:
        return False
    if not is_prime(q):
        return False 
    if not is_prime(r):
        return False
    if not fundamental_discriminant(D) == D:
        return False
    if D % 4 == 0: #check CM equation
        if not is_square(4*(t*t - 4*q)//D):
            return False
    if D % 4 == 1:
        if not is_square((t*t - 4*q)//D):
            return False
    if not (q+1-t) % r == 0: #check r | #E(F_q)
        return False
    if not power_mod(q,k,r) == 1: #check embedding degree is k
        return False
    return True


def filter_decorator(f):
    def helper(*args):
        q,t,r,k,D = f(*args)
        num_bits = _number_of_bits(r)
        while not is_suitable_curve(q,t,r,k,D, num_bits):
            q,t,r,k,D = f(*args)
            num_bits = _number_of_bits(r)
        return q,t,r,k,D
    return helper


def _number_of_bits(n):
    """
    Description:
        
        Returns the number of bits in the binary representation of n
    
    Input:
    
        n - integer
    
    Output:
    
        num_bits - number of bits
    
    """    
    if n == 0:
        return 1
    else:
        return floor(log(n).n()/log(2).n()) + 1


def is_suitable_curve(q,t,r,k,D, num_bits):
    """
    Description:
    
        User-defined method that filters the set of curves that are returned. By default checks if (q,t,r,k,D) is a valid curve
    
    Input:
    
        q - size of prime field
        t - trace of Frobenius
        r - size of prime order subgroup
        k - embedding degree
        D - (negative) fundamental discriminant
        num_bits - desired number of bits in r
    
    Output:
    
        bool - true iff (q,t,r,k,D) is suitable
    
    """
    return _number_of_bits(r) >= num_bits and is_valid_curve(q,t,r,k,D)


def is_suitable_q(q):
    """
    Description:
    
        User-defined method that filters the set of primes q that are returned. By default checks if q is prime
    
    Input:
    
        q - integer
    
    Output:
    
        bool - true iff q is suitable
    
    """
    return is_prime(q)


def is_suitable_r(r):
    """
    Description:
    
        User-defined method that filters the set of primes r that are returned. By default checks if r is prime
    
    Input:
    
        r - integer
    
    Output:
    
        bool - true iff r is suitable
    
    """
    return is_prime(r)
    
    
def print_curve(q,t,r,k,D):
    """
    Description:
    
        Prints the curve (q,t,r,k,D)
    
    Input:
    
        q - size of prime field
        t - trace of Frobenius
        r - size of prime order subgroup
        k - embedding degree
        D - (negative) fundamental discriminant
    
    Output:
    
        None
    
    """
    print(curve_to_string(q,t,k,r,D))
    
def curve_to_string(q,t,k,r,D):
    """
    Description:
    
        Returns a string representation of the curve (q,t,r,k,D)
    
    Input:
    
        q - size of prime field
        t - trace of Frobenius
        r - size of prime order subgroup
        k - embedding degree
        D - (negative) fundamental discriminant
    
    Output:
    
        s - string representation of (q,t,r,k,D)
    
    """
    if q == 0 or t == 0 or r == 0 or k == 0 or D == 0:
        return 'Failed to find an elliptic curve'
    else:
        return 'Elliptic curve over a field of size ' + str(q) + ' with trace ' + str(t) + ', a subgroup of order ' + str(r) + ' with embedding degree ' + str(k) + ', and fundamental discriminant ' + str(D)


def find_element_of_order(k,r):
    """
    Description:
    
        Finds a random element of order k in Z_r^*
    
    Input:
    
        k - integer such that r % k == 1
        r - prime
    
    Output:
    
        h - element of order k in Z_r^*
    
    """
    assert r % k == 1
    h = 0
    def order(h,k,p):
        bool = True
        g = h
        for i in range(1,k):
            bool = bool and g != 1
            g = Mod(g * h, p)
        bool = bool and g == 1
        return bool
    while not order(h,k,r): # expected number of iterations is k/euler_phi(k)
        h = power_mod(randint(2, r-1), (r-1)//k, r)
    return h


def method(r,k,D,max_trials=10000, g=0): 
    """
    Description:
        
        Run the Cocks-Pinch method to find an elliptic curve
    
    Input:
    
        r - prime
        k - embedding degree, r % k == 1
        D - (negative) fundamental discriminant where D is a square mod r
        max_trials - the number of integers q to test for primality in the CP method
        g - an element of order k in Z_r^*
    Output:
    
        (q,t) - tuple where q is a prime and t is chosen such that there exists
                an elliptic curve E over F_q with trace t, and r | q+1-t;
                if the algorithm fails to find (q,t), it will return (0,0)
    
    """
    assert test_promise(r,k,D), 'Invalid inputs'
    if g != 0:
        assert power_mod(g,k,r) == 1, 'Invalid inputs'
    else:
        g = find_element_of_order(k,r)
    D = Integer(D)
    t = Integer(g) + 1
    root_d = Integer(Mod(D, r).sqrt())
    u = Integer(Mod((t-2)*root_d.inverse_mod(r) ,r))
    q = 1
    j = Integer(0)
    i = Integer(0)
    count = 0
    while (count < max_trials):
           q = Integer( (t+i*r)**2 - D*(u + j*r)**2)
           if q % 4 ==0:
                q = q//4
                if is_suitable_q(q):
                    return (q, t+i*r)
                q = 1
           if random() < 0.5:
                j+=1
           else:
                i+=1
           count+=1
    return (0, 0) # no prime found, so end


@filter_decorator
def run(r,k,D,max_run_time=20):
    """
    Description:
    
        Runs the Cocks-Pinch method multiple times until a valid curve is found
    
    Input:
    
        r - prime
        k - embedding degree, r % k == 1
        D - (negative) fundamental discriminant where D is a square mod r
        max_run_time - maximum runtime of the function, in seconds
    Output:
        
        (q,t,r,k,D) - elliptic curve
    
    """
    assert test_promise(r,k,D), 'Invalid inputs'
    q = 0
    t = 0
    start = time.time()
    while q == 0 and time.time() - start < max_run_time: # run cocks pinch method until successful
        q,t = method(r,k,D)
    assert is_valid_curve(q,t,r,k,D), 'Invalid output'
    return q,t,r,k,D


def gen_params_from_bits(num_bits,k): 
    """
    Description:
    
        Generates a prime r with num_bits bits and a fundamental discriminant D to use as input to the Cocks-Pinch method
    
    Input:
        
        num_bits - number of bits in r
        k - embedding degree
        
    Output:
        
        r - prime such that r % k == 1 and r is num_bits bits long
        k - embedding degree
        D - (negative) fundamental discriminant where D is a square mod r
    
    """
    r = random_prime(2**num_bits, lbound=2**(num_bits-1))
    while not (r % k == 1 and is_suitable_r(r)):
        r = random_prime(2**num_bits, lbound=2**(num_bits-1))
    return gen_params_from_r(r,k)


def gen_params_from_r(r,k):
    """
    Description:
    
        Finds a fundamental discriminant D to use as input to the Cocks-Pinch method
    
    Input:
    
        r - prime such that r % k == 1
        k - embedding degree  
    
    Output:
        
        r - prime such that r % k == 1
        k - embedding degree
        D - (negative) fundamental discriminant where D is a square mod r
    
    """
    D = -Integer(Mod(int(random()*(1000)),r))
    i = 0
    while not kronecker(D,r) == 1: # expected number of iterations of the while loop is 2
        D = -Integer(Mod(int(random()*(1000)),r))
        i+=1
    D = fundamental_discriminant(D)
    if not (kronecker(D,r) == 1):
        return r, k, 0
    return r,k,D


def test_promise(r,k,D):
    """
    Description:
    
        Tests that r,k,D is a valid input to the Cocks-Pinch method
    
    Input:
    
        r - prime
        k - embedding degree    
        D - (negative) funadmental discriminant
    
    Output:
    
        bool - true iff (r,k,D) is a valid input to the Cocks-Pinch method
    
    """
    bool = (kronecker(D,r) == 1) # D is a square mod r
    bool = bool and ( (r-1) % k ==0) # k | r-1
    bool = bool and (D == fundamental_discriminant(D)) # check that D is a fundamental discriminant
    return bool


def make_curve(q,t,r,k,D,debug=False):
    """
    Description:
    
        Finds the curve equation for the elliptic curve (q,t,r,k,D) using the Complex Multiplication method
    
    Input:
    
        q - size of prime field
        t - trace of Frobenius
        r - size of prime order subgroup
        k - embedding degree
        D - (negative) fundamental discriminant
    
    Output:
    
        E - elliptic curve over F_q with trace t,
            a subgroup of order r with embedding degree k,
            and fundamental discriminant D
    
    """
    assert is_valid_curve(q,t,r,k,D), 'Invalid input. No curve exists.' # check inputs
    if debug:
        print('Tested input')
    poly = hilbert_class_polynomial(D) # compute hilbert class polynomial
    if debug:
        print('Computed Hilbert class polynomial')
    check = False
    j_inv = poly.any_root(GF(q)) # find j-invariant    
    orig_curve = EllipticCurve(GF(q), j=j_inv) # make a curve
    E = orig_curve
    check = test_curve(q,t,r,k,D,E) # see if this is the right curve
    twist = False
    if not check: # not the right curve, use quadratic twist
        E = E.quadratic_twist()
        check = test_curve(q,t,r,k,D,E)
        if check:
            twist = True
        else: # twist didnt work => j = 0 or 1728
            if j_inv == 0: # for j = 0, use sextic twists
                prim = primitive_root(q)
                i = 1
                while t != E.trace_of_frobenius() and i < 6:
                    E = orig_curve.sextic_twist(power_mod(prim,i,q))
                    i+=1
            elif j_inv == 1728: # for j = 1728, use quartic twists
                prim = primitive_root(q)
                i = 1
                while t != E.trace_of_frobenius() and i < 4:
                    E = orig_curve.quartic_twist(power_mod(prim,i,q))
                    i+=1
            else: # twist didnt work and j != 0, 1728. this should never happen
                raise Exception('Error. Quadratic twist failed to find the correct curve with j != 0, 1728')
            check = test_curve(q,t,r,k,D,E)
            twist = True
    if not check: # didnt find a curve. this should never happen, so write input to a file for debugging
        raise Exception('Error. Failed to find curve. Logging output to debug.txt')
        return False
    return E


def test_curve(q,t,r,k,D,E): 
    """
    Description:
    
       Tests that E is an elliptic curve over F_q with trace t, a subgroup of order r with embedding degree k, and fundamental discriminant D
    
    Input:
    
        q - size of prime field
        t - trace of Frobenius
        r - size of prime order subgroup
        k - embedding degree
        D - (negative) fundamental discriminant
    
    Output:
    
        bool - true iff E is an elliptic curve over F_q with trace t, a subgroup of order r with embedding degree k, and fundamental discriminant D
    
    """    
    bool = True
    bool = bool and (power_mod(q, k, r) == 1) #q^k -1 ==0 mod r
    bool = bool and (E.trace_of_frobenius() == t)
    bool = bool and (kronecker((t*t-4*q) * Integer(D).inverse_mod(q),q) == 1)
    bool = bool and (E.cardinality() == q+1-t)
    bool = bool and (E.cardinality() % r ==0)
    return bool


def get_curve_params(E, r, k):
    
    base_field = E.base_field()
    
    
    #find generator of G1
    
    G1 = E.random_point()
    while G1.order() % r != 0:
        G1 = E.random_point()
    G1 = G1 * int(G1.order() / r)

    if G1.order() != r:
        raise Exception("G1 Error!")
        return False
    
    #check embedding degree: we do only work for extensions of degree 4 and 6
    
    assert k == 4 or k == 6
    d = k/2
    
    #construct Extension field:
    #for k == 4 as a tower of degrees 2 : 2
    #for k == 6 as a tower of degrees 3 : 2
    # each extension is constructed via irreducible polynomial of the form X^n - \nu
    # where n is 2 or 3 and \nu is a quadratic (respectively cubic) nonresidue
    
    if k == 4:
        nu = base_field.random_element()
        while nu.is_square():
            nu = base_field.random_element()
        R.<x> = base_field[]
        twist_field = base_field.extension(x^2 - nu)
        assert twist_field.order() == base_field.order() ^ 2
    else:
        if base_field.order() % 3 != 1:
            raise Excepiion("No cubic non-residue in field!")
        
        gen = base_field.multiplicative_generator()
        x = randint(1, base_field.order() - 1)
        while x % 3 == 0:
            x = randint(1, base_field.order() - 1)
        nu = gen ^ x
        R.<x> = base_field[]
        twist_field = base_field.extension(x^3 - nu)
        assert twist_field.order() == base_field.order() ^ 3       
    first_nonresidue = nu
    
    nu = twist_field.random_element()
    while nu.is_square():
        nu = base_field.random_element()
    second_nonresidue = nu
    
    R.<y> = twist_field[]
    extension_field = twist_field.extension(x^2 - nu)
    assert extension_field.order() == twist_field.order() ^ 2
    
    twist = second_nonresidue
    
    #get the twisted curve
    #twisted short Weierstrass curve E'/twist_field : y^2 = x^3 + (a * twist^2) * x + (b * twist^3)
    
    twisted_curve = EllipticCurve(twist_field, [E.A * twist^2, E.b * twist^3])
    np = twisted_curve.cardinality()
    
    #get the generator of G2
    
     
    G1 = E.random_point()
    while G1.order() % r != 0:
        G1 = E.random_point()
    G1 = G1 * int(G1.order() / r)

    if G1.order() != r:
        raise Exception("G1 Error!")
        return False



k = 7 # embedding degree
num_bits = 100 # number of bits in size of prime order subgroup
r,k,D = gen_params_from_bits(num_bits,k)
q,t,r,k,D = run(r,k,D) # use CP method to solve for q and t
E = make_curve(q,t,r,k,D)


import time

    
def is_valid_curve(q,t,r,k,D): 
    """
    Description:
    
        Tests that (q,t,r,k,D) is a valid elliptic curve
    
    Input:
    
        q - size of prime field
        t - trace of Frobenius
        r - size of prime order subgroup
        k - embedding degree
        D - (negative) fundamental discriminant
    
    Output:
    
        bool - true iff there exists an elliptic curve over F_q with trace t, a subgroup of order r with embedding degree k, and fundamental discriminant D
    
    """
    if q == 0 or t == 0 or r == 0 or k == 0 or D == 0:
        return False
    if not is_prime(q):
        return False 
    if not is_prime(r):
        return False
    if not fundamental_discriminant(D) == D:
        return False
    if D % 4 == 0: #check CM equation
        if not is_square(4*(t*t - 4*q)//D):
            return False
    if D % 4 == 1:
        if not is_square((t*t - 4*q)//D):
            return False
    if not (q+1-t) % r == 0: #check r | #E(F_q)
        return False
    if not power_mod(q,k,r) == 1: #check embedding degree is k
        return False
    return True


def filter_decorator(f):
    def helper(*args):
        q,t,r,k,D = f(*args)
        num_bits = _number_of_bits(r)
        while not is_suitable_curve(q,t,r,k,D, num_bits):
            q,t,r,k,D = f(*args)
            num_bits = _number_of_bits(r)
        return q,t,r,k,D
    return helper


def _number_of_bits(n):
    """
    Description:
        
        Returns the number of bits in the binary representation of n
    
    Input:
    
        n - integer
    
    Output:
    
        num_bits - number of bits
    
    """    
    if n == 0:
        return 1
    else:
        return floor(log(n).n()/log(2).n()) + 1


def is_suitable_curve(q,t,r,k,D, num_bits):
    """
    Description:
    
        User-defined method that filters the set of curves that are returned. By default checks if (q,t,r,k,D) is a valid curve
    
    Input:
    
        q - size of prime field
        t - trace of Frobenius
        r - size of prime order subgroup
        k - embedding degree
        D - (negative) fundamental discriminant
        num_bits - desired number of bits in r
    
    Output:
    
        bool - true iff (q,t,r,k,D) is suitable
    
    """
    return _number_of_bits(r) >= num_bits and is_valid_curve(q,t,r,k,D)


def is_suitable_q(q):
    """
    Description:
    
        User-defined method that filters the set of primes q that are returned. By default checks if q is prime
    
    Input:
    
        q - integer
    
    Output:
    
        bool - true iff q is suitable
    
    """
    return is_prime(q)


def is_suitable_r(r):
    """
    Description:
    
        User-defined method that filters the set of primes r that are returned. By default checks if r is prime
    
    Input:
    
        r - integer
    
    Output:
    
        bool - true iff r is suitable
    
    """
    return is_prime(r)
    
    
def print_curve(q,t,r,k,D):
    """
    Description:
    
        Prints the curve (q,t,r,k,D)
    
    Input:
    
        q - size of prime field
        t - trace of Frobenius
        r - size of prime order subgroup
        k - embedding degree
        D - (negative) fundamental discriminant
    
    Output:
    
        None
    
    """
    print(curve_to_string(q,t,k,r,D))
    
def curve_to_string(q,t,k,r,D):
    """
    Description:
    
        Returns a string representation of the curve (q,t,r,k,D)
    
    Input:
    
        q - size of prime field
        t - trace of Frobenius
        r - size of prime order subgroup
        k - embedding degree
        D - (negative) fundamental discriminant
    
    Output:
    
        s - string representation of (q,t,r,k,D)
    
    """
    if q == 0 or t == 0 or r == 0 or k == 0 or D == 0:
        return 'Failed to find an elliptic curve'
    else:
        return 'Elliptic curve over a field of size ' + str(q) + ' with trace ' + str(t) + ', a subgroup of order ' + str(r) + ' with embedding degree ' + str(k) + ', and fundamental discriminant ' + str(D)


def find_element_of_order(k,r):
    """
    Description:
    
        Finds a random element of order k in Z_r^*
    
    Input:
    
        k - integer such that r % k == 1
        r - prime
    
    Output:
    
        h - element of order k in Z_r^*
    
    """
    assert r % k == 1
    h = 0
    def order(h,k,p):
        bool = True
        g = h
        for i in range(1,k):
            bool = bool and g != 1
            g = Mod(g * h, p)
        bool = bool and g == 1
        return bool
    while not order(h,k,r): # expected number of iterations is k/euler_phi(k)
        h = power_mod(randint(2, r-1), (r-1)//k, r)
    return h


def method(r,k,D,max_trials=10000, g=0): 
    """
    Description:
        
        Run the Cocks-Pinch method to find an elliptic curve
    
    Input:
    
        r - prime
        k - embedding degree, r % k == 1
        D - (negative) fundamental discriminant where D is a square mod r
        max_trials - the number of integers q to test for primality in the CP method
        g - an element of order k in Z_r^*
    Output:
    
        (q,t) - tuple where q is a prime and t is chosen such that there exists
                an elliptic curve E over F_q with trace t, and r | q+1-t;
                if the algorithm fails to find (q,t), it will return (0,0)
    
    """
    assert test_promise(r,k,D), 'Invalid inputs'
    if g != 0:
        assert power_mod(g,k,r) == 1, 'Invalid inputs'
    else:
        g = find_element_of_order(k,r)
    D = Integer(D)
    t = Integer(g) + 1
    root_d = Integer(Mod(D, r).sqrt())
    u = Integer(Mod((t-2)*root_d.inverse_mod(r) ,r))
    q = 1
    j = Integer(0)
    i = Integer(0)
    count = 0
    while (count < max_trials):
           q = Integer( (t+i*r)**2 - D*(u + j*r)**2)
           if q % 4 ==0:
                q = q//4
                if is_suitable_q(q):
                    return (q, t+i*r)
                q = 1
           if random() < 0.5:
                j+=1
           else:
                i+=1
           count+=1
    return (0, 0) # no prime found, so end


@filter_decorator
def run(r,k,D,max_run_time=20):
    """
    Description:
    
        Runs the Cocks-Pinch method multiple times until a valid curve is found
    
    Input:
    
        r - prime
        k - embedding degree, r % k == 1
        D - (negative) fundamental discriminant where D is a square mod r
        max_run_time - maximum runtime of the function, in seconds
    Output:
        
        (q,t,r,k,D) - elliptic curve
    
    """
    assert test_promise(r,k,D), 'Invalid inputs'
    q = 0
    t = 0
    start = time.time()
    while q == 0 and time.time() - start < max_run_time: # run cocks pinch method until successful
        q,t = method(r,k,D)
    assert is_valid_curve(q,t,r,k,D), 'Invalid output'
    return q,t,r,k,D


def gen_params_from_bits(num_bits,k): 
    """
    Description:
    
        Generates a prime r with num_bits bits and a fundamental discriminant D to use as input to the Cocks-Pinch method
    
    Input:
        
        num_bits - number of bits in r
        k - embedding degree
        
    Output:
        
        r - prime such that r % k == 1 and r is num_bits bits long
        k - embedding degree
        D - (negative) fundamental discriminant where D is a square mod r
    
    """
    r = random_prime(2**num_bits, lbound=2**(num_bits-1))
    while not (r % k == 1 and is_suitable_r(r)):
        r = random_prime(2**num_bits, lbound=2**(num_bits-1))
    return gen_params_from_r(r,k)


def gen_params_from_r(r,k):
    """
    Description:
    
        Finds a fundamental discriminant D to use as input to the Cocks-Pinch method
    
    Input:
    
        r - prime such that r % k == 1
        k - embedding degree  
    
    Output:
        
        r - prime such that r % k == 1
        k - embedding degree
        D - (negative) fundamental discriminant where D is a square mod r
    
    """
    D = -Integer(Mod(int(random()*(1000)),r))
    i = 0
    while not kronecker(D,r) == 1: # expected number of iterations of the while loop is 2
        D = -Integer(Mod(int(random()*(1000)),r))
        i+=1
    D = fundamental_discriminant(D)
    if not (kronecker(D,r) == 1):
        return r, k, 0
    return r,k,D


def test_promise(r,k,D):
    """
    Description:
    
        Tests that r,k,D is a valid input to the Cocks-Pinch method
    
    Input:
    
        r - prime
        k - embedding degree    
        D - (negative) funadmental discriminant
    
    Output:
    
        bool - true iff (r,k,D) is a valid input to the Cocks-Pinch method
    
    """
    bool = (kronecker(D,r) == 1) # D is a square mod r
    bool = bool and ( (r-1) % k ==0) # k | r-1
    bool = bool and (D == fundamental_discriminant(D)) # check that D is a fundamental discriminant
    return bool


def make_curve(q,t,r,k,D,debug=False):
    """
    Description:
    
        Finds the curve equation for the elliptic curve (q,t,r,k,D) using the Complex Multiplication method
    
    Input:
    
        q - size of prime field
        t - trace of Frobenius
        r - size of prime order subgroup
        k - embedding degree
        D - (negative) fundamental discriminant
    
    Output:
    
        E - elliptic curve over F_q with trace t,
            a subgroup of order r with embedding degree k,
            and fundamental discriminant D
    
    """
    assert is_valid_curve(q,t,r,k,D), 'Invalid input. No curve exists.' # check inputs
    if debug:
        print('Tested input')
    print "hilbert class poly"
    poly = hilbert_class_polynomial(D) # compute hilbert class polynomial
    print "hilbert poly found"
    if debug:
        print('Computed Hilbert class polynomial')
    check = False
    print "j_inv"
    j_inv = poly.any_root(GF(q)) # find j-invariant 
    print "j_inv found"
    orig_curve = EllipticCurve(GF(q), j=j_inv) # make a curve
    E = orig_curve
    check = test_curve(q,t,r,k,D,E) # see if this is the right curve
    twist = False
    if not check: # not the right curve, use quadratic twist
        E = E.quadratic_twist()
        check = test_curve(q,t,r,k,D,E)
        if check:
            twist = True
        else: # twist didnt work => j = 0 or 1728
            if j_inv == 0: # for j = 0, use sextic twists
                prim = primitive_root(q)
                i = 1
                while t != E.trace_of_frobenius() and i < 6:
                    E = orig_curve.sextic_twist(power_mod(prim,i,q))
                    i+=1
            elif j_inv == 1728: # for j = 1728, use quartic twists
                prim = primitive_root(q)
                i = 1
                while t != E.trace_of_frobenius() and i < 4:
                    E = orig_curve.quartic_twist(power_mod(prim,i,q))
                    i+=1
            else: # twist didnt work and j != 0, 1728. this should never happen
                raise Exception('Error. Quadratic twist failed to find the correct curve with j != 0, 1728')
            check = test_curve(q,t,r,k,D,E)
            twist = True
    if not check: # didnt find a curve. this should never happen, so write input to a file for debugging
        raise Exception('Error. Failed to find curve. Logging output to debug.txt')
        return False
    return E


def test_curve(q,t,r,k,D,E): 
    """
    Description:
    
       Tests that E is an elliptic curve over F_q with trace t, a subgroup of order r with embedding degree k, and fundamental discriminant D
    
    Input:
    
        q - size of prime field
        t - trace of Frobenius
        r - size of prime order subgroup
        k - embedding degree
        D - (negative) fundamental discriminant
    
    Output:
    
        bool - true iff E is an elliptic curve over F_q with trace t, a subgroup of order r with embedding degree k, and fundamental discriminant D
    
    """    
    bool = True
    bool = bool and (power_mod(q, k, r) == 1) #q^k -1 ==0 mod r
    bool = bool and (E.trace_of_frobenius() == t)
    bool = bool and (kronecker((t*t-4*q) * Integer(D).inverse_mod(q),q) == 1)
    bool = bool and (E.cardinality() == q+1-t)
    bool = bool and (E.cardinality() % r ==0)
    return bool


def get_final_exp_params(x, p):
    c = int(x/p)
    d = x - p*c
    if d > (p/2):
        d -= p
        c += 1
    if c >= p:
        return False
    else:
        return (c, d)


def get_curve_params(E, r, k, max_run_time=20):
    
    base_field = E.base_field()
    
    print "first point count"
    point_count = E.cardinality()
    print "end first point count"
    assert point_count % r == 0
    h = int(point_count / r)
    
    #find generator of G1
    print "finding point"
    G1 = h * E.random_point()
    while G1 == E([0,1,0]):
        G1 = h * E.random_point()
    print "finding point end"

    #check embedding degree: we do only work for extensions of degree 4 and 6
    
    assert k == 4 or k == 6
    d = k/2
    
    #construct Extension field:
    #for k == 4 as a tower of degrees 2 : 2
    #for k == 6 as a tower of degrees 3 : 2
    # each extension is constructed via irreducible polynomial of the form X^n - \nu
    # where n is 2 or 3 and \nu is a quadratic (respectively cubic) nonresidue
    
    print "finding extension generators"
    
    if k == 4:
        flag = False
        start = time.time()
        while not flag and time.time() - start < max_run_time:
            nu = base_field.random_element()
            while nu.is_square():
                nu = base_field.random_element()
            R.<x> = base_field[]
            twist_field.<u> = base_field.extension(x^2 - nu)
            if not u.is_square():
                flag = true
                non_residue = nu
                twist = u
        if not flag:
            raise Exception("No suitable degree tower")
    else:
        if base_field.order() % 3 != 1:
            raise Exception("No cubic non-residue in field!")
        
        gen = base_field.multiplicative_generator()
        flag = False
        start = time.time()
        while not flag and time.time() - start < max_run_time:
            x = randint(1, base_field.order() - 1)
            while x % 3 == 0:
                x = randint(1, base_field.order() - 1)
            nu = gen ^ x
            R.<x> = base_field[]
            twist_field.<u> = base_field.extension(x^3 - nu)
            if not u.is_square():
                flag = true
                non_residue = nu
                twist = u
        if not flag:
            raise Exception("No suitable degree tower")
    
    #get the twisted curve
    #twisted short Weierstrass curve E'/twist_field : y^2 = x^3 + (a * twist^2) * x + (b * twist^3)
    print "finding extension generators end."
    twisted_curve = EllipticCurve(twist_field, [E.a4() * twist^2, E.a6() * twist^3])
    
    print "second point count"
    n = twisted_curve.cardinality()
    print "end second point count"
    assert n % r == 0
    h = int(n / r)
    
    #get the generator of G2
    
    print "second point gen"     
    G2 = h * twisted_curve.random_point()
    while G2 == twisted_curve([0, 1, 0]):
        G2 = h * twisted_curve.random_point()
    print "end second point gen"
    
    #additional parameters: length of Miller loop and data for final exponentiation
    p = base_field.order()
    T = p - point_count
    
    if k == 4:
        #hard part of final exponentiation: (p^2+1)/r = w_1 * p + w_0
        assert (p^2 + 1) % r == 0
        a = int((p^2+1)/r)
    else:
        # (p^2-p+1)/r = w_1 * p + w_0
        assert (p^2-p+1) % r == 0
        a = int((p^2-p+1)/r)
    
    exp_params = get_final_exp_params(a, p)
    if not exp_params:
            raise Exception("Incorrect final exp decomposition")
    w_1, w_0 = exp_params
    
    return p, non_residue, twisted_curve, G1, G2, T, w_1, w_0
    
    
def print_all_data(p, r, k, nonresidue, curve, twisted_curve, G1, G2, T, w_1, w_0, stream):
    data = {}
    data["base_field"] = str(p)
    data["base field bitlength"] = str(p.nbits())
    data["group order"] = str(r)
    data["group order bitlength"] = str(r.nbits())
    data["embedding degree"] = str(k)
    data["extension non-residue"] = str(nonresidue)
    data["curve A"] = str(curve.a4())
    data["curve B"] = str(curve.a6())
    data["twisted curve A"] = str(twisted_curve.a4())
    data["twisted curve B"] = str(twisted_curve.a6())
    data["generator of first source group"] = str(G1)
    data["generator of second source group"] = str(G2)
    data["length of Miller loop"] = str(T)
    data["w1 coef for final exp"] = str(w_1)
    data["w0 coed for final exp"] = str(w_0)
    
    stream.write(json.dumps(data, indent=4, separators=(',', ': ')))
    stream.write("\n")
    
    
def generate_test_curves(filename):
    f = open(filename,'w')
    for k in (4, 6):
        for num_bits in xrange(325, 400, 25):
            found = False
            print "..........."
            while not found:
                try:
                    print "1"
                    r,k,D = gen_params_from_bits(num_bits,k)
                    print "2"
                    q,t,r,k,D = run(r,k,D) # use CP method to solve for q and t
                    print "3"
                    E = make_curve(q,t,r,k,D)
                    print "4"
                    p, non_residue, twisted_curve, G1, G2, T, w_1, w_0 = get_curve_params(E, r, k)
                    print "5"
                    print_all_data(p, r, k, non_residue, E, twisted_curve, G1, G2, T, w_1, w_0, f)
                    print "6"
                    found = True
                except Exception as e:
                    print str(e)
    f.close()


def test_func():
    k = 6 # embedding degree
    num_bits = 10 # number of bits in size of prime order subgroup
    r,k,D = gen_params_from_bits(num_bits,k)
    q,t,r,k,D = run(r,k,D) # use CP method to solve for q and t
    E = make_curve(q,t,r,k,D)
    p, non_residue, twisted_curve, G1, G2, T, w_1, w_0 = get_curve_params(E, r, k)
    print_all_data(p, r, k, non_residue, E, twisted_curve, G1, G2, T, w_1, w_0, sys.stdout)

    
#generate_test_curves("ololo.txt")
generate_test_curves("/home/k/curves2.txt")

                 

