p = 47
R = Zp(p, prec = 3, type = 'fixed-mod', print_mode = 'terse')

def find_coefs():
    n = p^2
    roots = []
    for x in xrange(n):
        if (x^3 - 2*x^2+4*x-4) % n == 0 and 3*x^2 - 4*x + 4 % p !=0:
            if all([y - x % n for y in roots]):
                roots.append(x)
    return roots


gamma, beta, alpha =  (R(x) for x in find_coefs())

MR = MatrixSpace(R,3,3)
main_mat = MR([1, 1, 1, alpha, beta, gamma, alpha^2, beta^2, gamma^2])
d = main_mat.det()
A = (gamma - beta) / d
B = (alpha - gamma) / d
C = (beta - alpha) / d
print A, B, C

a = alpha^46 - 1
b = beta^46 - 1
c = gamma^46 - 1

possible_r = []
for r in xrange(p-1):
    x = A*alpha^r + B*beta^r+C*gamma^r
    if not x[0]:
        possible_r.append(r)
        
print possible_r

s = 1
for r in (0, 1, 4, 6, 13):
    x = A*alpha^r*(1+a)^s + B*beta^r*(1+b)^s + C*gamma^r*(1+c)^s
    print x[0], x[1], x[2]
    
r = 6 