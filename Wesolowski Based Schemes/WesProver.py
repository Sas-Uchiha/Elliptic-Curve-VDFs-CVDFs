'''Input: Command Line Inputs: T, the delay parameter followed by x: unhashed x-coordinate of the starting point '''
'''Dependencies: EC.py'''
import EC as elliptic_curve
import sys
from hashlib import sha256
import time
import math

def setup(T):
    a,b,p = elliptic_curve.give_parameters()
    ec = elliptic_curve.EllipticCurve(a,b,p)
    return ec #We are implicitly assuming that we are using SHA256.

def evaluation(pp, starting_point):
    ec = pp[0]
    T = pp[1]
    #evaluation start here
#    start = time.time()
    y = ec.multiply_point_by_2T(T, starting_point, 1)
#    end = time.time()
#    print(f"evaluation time: {round((end - start) * 1000,2)}")
    #evaluation stop here
    return y
    
def highestPowerof2(n):
    res = 0;
    for i in range(n, 0, -1):
        if ((i & (i - 1)) == 0):
            res = i;
            break;
    return res;

def optimised_power(T, g, curve, l, sequence):
    k = highestPowerof2(int(int(math.log(T)//math.log(2)) / 2))
    gamma = 1
    k1 = k//2
    k0 = k - k1
    x = (-1,-1)
    power = 2**k
    power0 = 2**k0
#    C = [curve.multiply_point_by_2T(i*k*gamma, g) for i in range(0, T//(k * gamma) + 1)]
#    if math.ceil(math.log(T)/math.log(2)) == math.floor(math.log(T)/math.log(2)):
    C = [sequence[i * k * gamma] for i in range(0, math.ceil(T//(k * gamma)) + 1)]
#    else:
#        C = [sequence[i * k * gamma] for i in range(0, math.ceil(T//(k * gamma)))]
#        point = C[-1]
#        C.append(curve.multiply_point_by_k(2**(k*gamma), point))
#    print(C1 == C)
    for j in range(gamma - 1, -1, -1):
        x = curve.multiply_point_by_2T(k, x)
        y = [(-1,-1) for i in range(power)]
        for i in range(0, math.ceil(T/(k * gamma))):
            b = get_block(k, l, T, i*gamma + j)
            y[b] = curve.addition(y[b], C[i])
        for b1 in range(0, 2**k1):
            z = (-1,-1)
            for b0 in range(0, power0):
                z = curve.addition(z, y[b1*power0 + b0])
            x = curve.addition(x, curve.multiply_point_by_k(b1 * power0, z))
        for b0 in range(0, power0):
            z = (-1,-1)
            for b1 in range(0, 2**k1):
                z = curve.addition(z, y[b1*power0 + b0])
            x = curve.addition(x, curve.multiply_point_by_k(b0, z))
    return x

def get_block(k, l, t, i):
    return (2**k) * pow(2, t - k * (i + 1), l) // l

def verify(T, curve, g, h, sequence):
    if not curve.check_point(g) or not curve.check_point(h):
        return False
        
    #proof time prover start here
    
    l = elliptic_curve.generate_prime_number_below_bits(133)
#     print("Size = " + str(sys.getsizeof(l)))
#    start = time.time()
    pi = optimised_power(T, g, curve, l, sequence)
#    end = time.time()
#    print(f"proof time: {round((end - start) * 1000,2)}")
    #proof time prover stop here
    #proof time verifier start here
    r = pow(2, T, l)
#     start = time.time()
    final_point = curve.addition(curve.multiply_point_by_k(l, pi), curve.multiply_point_by_k(r, g))
    if h == final_point:
        b = True
    else:
        b = False
#     end = time.time()
#     print(f"verification time: {round((end - start) * 1000,2)}")
    #proof time verifier stop here
    return b

if __name__ == "__main__":
    T = int(sys.argv[1])
    x = sys.argv[2]
    curve = setup(T)
    print(curve)
    g = int(sha256(x.encode('utf-8')).hexdigest(),16)
    starting_point = curve.curve_point(g)
    coordinates = evaluation((curve, T), starting_point)
    h = coordinates[0]
    sequence = coordinates[1]
    if not h[0] == h[1] == -1:
        print(verify(T, curve, starting_point, h, sequence))
