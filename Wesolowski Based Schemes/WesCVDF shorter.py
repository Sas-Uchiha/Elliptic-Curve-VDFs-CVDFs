'''Input: Command Line Inputs: T, the delay parameter followed by t: time upto which the cVDF needs to be calculated'''
'''Dependencies: EC.py, WesProver.py, WesUVDF.py''' 

import WesProver as wes
import EC as elliptic_curve
import sys
from hashlib import sha256
import math as math
import WesUVDF as uvdf
import time

k = 8 #change here if required
b = 2
power = k**b

def setup(T):
    return uvdf.setup(T)
    
def generate():
    return uvdf.generate()

def eval(curve, x):
    to_compute = t//power
    proof = ((-1, -1), 0, 0)
    sequence = [x]
    proof_time = 0
    for i in range(1, to_compute + 1):
        (x, seq) = curve.multiply_point_by_2T(power, x, 1)
        sequence.extend(seq[1:])
        l = elliptic_curve.generate_prime_number_below_bits(133)
        start = time.time()
        proof = (wes.optimised_power(i * power, x, curve, l, sequence), l, i * power)
        end = time.time()
        proof_time += end -start
    x = curve.multiply_point_by_2T(t % power, x)
    return (x, proof, proof_time)
    
def verify(curve, h, g, proof, t):
    if not curve.check_point(g) or not curve.check_point(h):
        return False
    if proof[0] == (-1, -1):
        return False
    delay = proof[2]
    l = proof[1]
    r = pow(2, delay, l)
    proof_till_multiple_of_kb = curve.addition(curve.multiply_point_by_k(l, proof[0]), curve.multiply_point_by_k(r, g))
    return curve.multiply_point_by_2T(t % power, proof_till_multiple_of_kb) == h

    
if __name__ == "__main__":
    T = int(sys.argv[1])
    t = int(sys.argv[2])
    curve = setup(T)
    print(curve)
    x0 = generate()
    while math.gcd(x0 - 1, curve.curve_params_prime) != 1:
        x0 = generate()
    x = curve.curve_point(x0)
#     start = time.time()
    (h, proof, time1) = eval(curve, x)
#     end = time.time()
#     print(f"Proof Time: {round(time1 * 1000,2)}")
#     print(f"Eval Time: {round((end - start - time1) * 1000,2)}")
#     print(f"Proof Size: {sys.getsizeof(proof)}")
    print(h)
    if not h == (-1, -1):
        start = time.time()
        b = verify(curve, h, x, proof, t)
        end = time.time()
#         print(f"Verify Time: {round((end - start) * 1000,2)}")
        print(b)
    
    
    
