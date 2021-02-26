'''Input: Command Line Inputs: T, the delay parameter followed by t: time upto which the cVDF needs to be calculated'''
'''Dependencies: EC.py, Pietr.py, PietrUVDF.py'''
import Pietr as pietr
import EC as elliptic_curve
import sys
from hashlib import sha256
import math as math
import PietrUVDF as uvdf
import time

k = 4 #change here if required
b = 1
T = 0
n = 0
class Node:
    def __init__(self, initial_point, final_point, xprime, yprime, level):
        self.initial_point = initial_point
        self.final_point = final_point
        self.xprime = xprime
        self.yprime = yprime
        self.children = []
        self.level = level
    
    def add_child(self, node):
        self.children.append(node)

def setup(T):
    return uvdf.setup(T)
    
def generate():
    return uvdf.generate()

def eval(curve, x, depth):
    if depth < b:
        return base_case(curve, x, depth)
    node = Node((-1,-1), (-1, -1), (-1, -1), (-1, -1), 0)
    node.children.append(eval(curve, x, depth -1))
    for i in range(1, k):
        child_node = eval(curve, node.children[-1].final_point, depth - 1)
        node.children.append(child_node)
    node.initial_point = node.children[0].initial_point
    node.final_point = node.children[-1].final_point
    node.level = node.children[0].level + 1
    node.xprime, node.yprime = FSHProver(curve, node)
    return node
    

def base_case(curve, x, depth):
    T1 = T//(k**(n - depth - 1))
#    print(n-depth -1)
    starting_point = x
    final_point = pietr.evaluation((curve, T1), starting_point)
    xprime = (-1, -1)
    yprime = (-1, -1)
    level = depth
    return Node(starting_point, final_point, xprime, yprime, level)
    

def FSHProver(curve, node):
    starting_points = [child.initial_point for child in node.children]
    final_points = [child.final_point for child in node.children]
    (xprime, yprime) = combine(curve, starting_points, final_points)
    return xprime, yprime
    
def combine(curve, starting_points, final_points):
    alpha = [elliptic_curve.generate_prime_number_below_bits(128) for i in range(k)]
    xprime = (-1,-1)
    yprime = (-1,-1)
    for i in range(k):
        xprime = curve.addition(xprime, curve.multiply_point_by_k(alpha[i], starting_points[i]))
    for i in range (k):
        yprime = curve.addition(yprime, curve.multiply_point_by_k(alpha[i], final_points[i]))
    return (xprime, yprime)
    
def verify(curve, node):
    if node.xprime == node.yprime == (-1, -1):
        return pietr.verify(T//(k**(n-node.level - 1)), curve, node.initial_point, node.final_point)
    status = True
    for child in node.children:
        status = status and verify(curve, child)
        if not status:
            break
    if not status:
        return status
    else:
        return FSHVerifier(T//(k**(n-node.level)),curve, node.xprime, node.yprime)
        
    
def FSHVerifier(depth,curve, xprime, yprime):
    return pietr.verify(depth, curve, xprime, yprime)

def numberToBase(n, b):
    if n == 0:
        return [0]
    digits = []
    while n:
        digits.append(int(n % b))
        n //= b
    return digits[::-1]

    
if __name__ == "__main__":
    T = int(sys.argv[1])
    t = int(sys.argv[2])
    curve = setup(T)
    print(curve)
    x0 = generate()
    while math.gcd(x0 - 1, curve.curve_params_prime) != 1:
        x0 = generate()
    time_in_base_k = numberToBase(t, k)
    print(time_in_base_k)
    tree = []
    n = len(time_in_base_k)
    x = curve.curve_point(x0)
    if x == (-1,-1):
        sys.exit()
#     start = time.time()
    for i in range(n - 1, -1, -1):
        for j in range(time_in_base_k[i]):
            node = eval(curve, x, n - i- 1)
            x = node.final_point
            tree.append(node)
    
#     end = time.time()
#     print(f"evaluation: {round((end-start) * 1000,2)}")
    status = True
    counter = 0
#     start = time.time()
    for node in tree:
        status = status & verify(curve, node)
        counter += k ** node.level
#     end = time.time()
#     print(tree[0].level)
#     print(f"verification: {round((end-start) * 1000,2)}")
    print(status and counter == t)
#     print(pietr.proof_time)
        
    
    
    
