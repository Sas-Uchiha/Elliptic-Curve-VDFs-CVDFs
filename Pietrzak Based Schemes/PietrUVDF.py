'''Input: Command Line Inputs: T, the delay parameter'''
'''Dependencies: EC.py, Pietr.py'''
import Pietr as pietr
import EC as elliptic_curve
import sys
from hashlib import sha256
import math as math

k = 4 #change here if required
b = 1

def setup(T):
    return pietr.setup(T)
    
def generate():
    return elliptic_curve.give_prime(128) + 1

def eval(curve, x, T):
    if math.gcd(x - 1, curve.curve_params_prime) != 1:
        return (-1, -1)
    starting_point = curve.curve_point(x)
    y = pietr.evaluation((curve, T), starting_point)
    if T < k ** b:
        return (y, -1)
    return (y, FSHProver(curve, starting_point, T))

def FSHProver(curve, x, T):
    (xprime, yprime) = combine(curve, x, T)
    #change below if required
    return (xprime, yprime)
    
def combine(curve, starting_point, T):
    seg = [starting_point]
    for i in range(1, k + 1):
        seg.append(curve.multiply_point_by_2T(i*T//k, starting_point))
    alpha = [elliptic_curve.generate_prime_number_below_bits(128) for i in range(k)]
    xprime = (-1,-1)
    yprime = (-1,-1)
    for i in range(1, k + 1):
        xprime = curve.addition(xprime, curve.multiply_point_by_k(alpha[i-1], seg[i-1]))
    for i in range (1, k + 1):
        yprime = curve.addition(yprime, curve.multiply_point_by_k(alpha[i-1], seg[i]))
    return (xprime, yprime)
    
def verify(curve, x, y, T, xprime, yprime):
    if not curve.check_point(x) or not curve.check_point(y):
        return False
    return pietr.verify(T//k, curve, xprime, yprime)
    
    
if __name__ == "__main__":
    T = int(sys.argv[1])
    x = generate()
    curve = setup(T)
    print(curve)
    starting_point = curve.curve_point(x)
    h = eval(curve, x, T)
    print(h)
    if type(h[1]) == int:
        print(pietr.verify(T, curve, starting_point, h[0]))
    else:
        print(verify(curve, starting_point, h[0], T, h[1][0], h[1][1]))
