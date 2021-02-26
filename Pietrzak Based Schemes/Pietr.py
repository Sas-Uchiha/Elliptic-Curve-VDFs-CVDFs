'''Input: Command Line Inputs: T, the delay parameter followed by x: unhashed x-coordinate of the starting point '''
'''Dependencies: EC.py'''

import EC as elliptic_curve
import sys
from hashlib import sha256
from random import randint
import time

proof_time = 0
proof_size = 0

def setup(T):
    a,b,p = elliptic_curve.give_parameters()
    ec = elliptic_curve.EllipticCurve(a,b,p)
    return ec #We are implicitly assuming that we are using SHA256.

def evaluation(pp, starting_point):
    ec = pp[0]
    T = pp[1]
    #evaluation start here
#    start = time.time()
    y = ec.multiply_point_by_2T(T, starting_point)
#    end = time.time()
#    print(f"evaluation time: {round((end - start)*1000,2)}")
    #evaluation stop here
    return y
    

def verify(T, curve, g, h):
#     global proof_time
    if T == 1:
        return h == curve.addition(g, g)
    #proof time analysis start here/proof size
#     start = time.time()
    v = curve.multiply_point_by_2T((T + 1)//2, g)
#     end = time.time()
#     proof_time += round((end - start)*1000,2)
#     proof_size += sys.getsizeof(v)
#    print(sys.getsizeof(v))
    #proof time analysis stop here/proof size
    if not curve.check_point(v):
        return False
    
    r = int(''.join(map(str,[randint(0, 1) for _ in range(0, 129)])),2)
    g1 = curve.addition(curve.multiply_point_by_k(r,g),v)
    h1 = curve.addition(h, curve.multiply_point_by_k(r,v))
    return verify( (T + 1) // 2, curve, g1, h1)

if __name__ == "__main__":
    T = int(sys.argv[1])
    x = sys.argv[2]
    curve = setup(T)
    print(curve)
    g = int(sha256(x.encode('utf-8')).hexdigest(),16)
    starting_point = curve.curve_point(g)
    h = evaluation((curve, T), starting_point)
    print(h)
#     start = time.time()
    #verification start here
    verified = verify(T, curve, starting_point, h)
#     end = time.time()
#     print(f"verification time: {round((end - start)*1000,2)}")
    #verification stop here subtract proof time.
    print(verified)
