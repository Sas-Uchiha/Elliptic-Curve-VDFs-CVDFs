'''No Dependencies'''
import hashlib
from random import SystemRandom, randint

class EllipticCurve:

    def __init__(self, a, b, p):
        self.p_length = 256
        self.curve_params_a = a
        self.curve_params_b = b
        self.curve_params_prime = p
    """

    point (-1, -1) is being treated as identity as for any curve such that 4a^3 + 27b^2 not = 0, the point (-1, -1) will never
     lie on the curve and hence is taken as the identity for convenience

    """

    """

    Calculate modular inverse of a number a with respect to a prime number :

    Input : n = number whose modular inverse is required
            p = prime number

    Output : inv = modular inverse of a

    """


    def mod_inv(self, n, p):
        return pow(n, p - 2, p)


    """

    Returns a^(p-1 // 2) mod p

    """


    def legendre(self, a, p):
        return pow(a, (p - 1) // 2, p)


    """

    Finds the Modular Square Root using Tonelli-Shank's Algorithm

    Input : n = number whose modular square root is to be found
            p = prime number

    Output : r = modular square root of n

    """


    def tonelli_shank(self, n, p):
        if self.legendre(n, p) == p - 1:
            return -1

        q = p - 1
        s = 0
        while q % 2 == 0:
            q //= 2
            s += 1

        if s == 1:
            return pow(n, (p + 1) // 4, p)

        z = 2

        while not self.legendre(z, p) == (p - 1):
            z += 1

        c = pow(z, q, p)
        r = pow(n, (q + 1) // 2, p)
        t = pow(n, q, p)
        m = s
        while True:

            if t == 0:
                return 0

            if t % p == 1:
                return r

            x = pow(t, 2, p)
            i = 1
            while not x == 1:
                x = pow(x, 2, p)
                i += 1

            if i == m:
                return -1

            b = pow(c, 1 << (m - i - 1), p)
            r = (r * b) % p
            c = (b * b) % p
            t = (t * c) % p
            m = i

    """

    Returns the order of a point in the curve

    Input : x: coordinates of a point
            p: prime modulus

    """

    def order(self, x):
        p = self.curve_params_prime
        count = 0
        k = x
        while k[0] != -1 and k[1] != -1:
            count += 1
            k = self.multiply_point_by_k(count, x)
            if count > 1e25:
                return 1

        return count



    """

    Returns the coordinates of an elliptic curve

    Input : x = x - coordinate
            a = coefficient of x
            b = constant
            p = prime number

    Output : (x, y) = coordinates of the point on the elliptic curve

    """


    def curve_point(self, x):
        p = self.curve_params_prime
        a = self.curve_params_a
        b = self.curve_params_b
        y = (x ** 3) + a * x + b
        y = self.tonelli_shank(y, p)
        if y == -1:
            return -1, -1
        return x, y

    def check_point(self, point):
        if point == (-1,-1):
            return True
        x = point[0]
        y = point[1]
        ysq = pow(y,2,self.curve_params_prime)
        rhs = (pow(x,3,self.curve_params_prime) + self.curve_params_b + self.curve_params_a * pow(x,1,self.curve_params_prime))% self.curve_params_prime
        return ysq == rhs
    """

    Performs Binary Operation of Addition on two given points

    Input : A and B = Points on the Elliptic Curve
            p = prime number

    Output : C = Result of A + B

    """


    def addition(self, x, y):
        p = self.curve_params_prime
        a = self.curve_params_a
        (x1, y1) = x
        (x2, y2) = y

        if x1 == y1 == -1:
            return y

        elif x2 == y2 == -1:
            return x

        if x1 == x2 and y1 != y2:
            return -1, -1

        if x != y:
            m = ((y2 - y1) * self.mod_inv(x2 - x1, p))

        else:
            m = (3 * (x1 ** 2) + a) * self.mod_inv(2 * y1, p)

        m %= p
        x3 = m ** 2 - x1 - x2
        x3 %= p
        y3 = (m * (x1 - x3) - y1) % p
        c = (x3, y3)

        return c



    def multiply_point_by_2T_helper(self, T, original_point):
        prime = self.curve_params_prime
        a = [original_point]

        binary_k = "1" + "0" * T
        length = T + 1

        for i in range(1, length):
            a.append(self.addition(a[i - 1], a[i - 1]))

        key = (-1, -1)

        for i in range(length):
            if binary_k[length - 1 - i] == '1':
                key = self.addition(key, a[i])

        return (key, a)
    
    def multiply_point_by_2T(self, T, original_point, params = 0):
        if params == 0:
            return self.multiply_point_by_2T_helper(T, original_point)[0]
        else:
            return self.multiply_point_by_2T_helper(T, original_point)

    def multiply_point_by_k(self, k, original_point):
        prime = self.curve_params_prime
        a = [original_point]

        binary_k = bin(k)[2:]
        length = len(binary_k)

        for i in range(1, length):
            a.append(self.addition(a[i - 1], a[i - 1]))

        key = (-1, -1)

        for i in range(length):
            if binary_k[length - 1 - i] == '1':
                key = self.addition(key, a[i])

        return key

    def __str__(self):
        return f"a = {self.curve_params_a}, b = {self.curve_params_b}, prime = {self.curve_params_prime}"

"""

generate a number of length

"""


def generate_number(k):
    random = SystemRandom()
    p = random.getrandbits(k)
    binary = bin(p)
    new_binary = binary[0:len(binary) - 1] + "1"
    return int(new_binary, 2)

def generate_prime_number_below_bits(x):
    prime = int(''.join(map(str,[randint(0,1) for _ in range(x)])),2)
    while not check_prime(prime):
        prime = int(''.join(map(str,[randint(0,1) for _ in range(x)])),2)

    return prime

"""

Performs Miller Rabin Test

Input : n = number, r

"""


def miller_rabin_test(r, n):
    a = 2 + randint(2, n - 2)

    p = pow(a, r, n)

    if p == 1 or p == n - 1:
        return True

    while r != (n - 1):

        p = pow(p, 2, n)
        r *= 2

        if p == 1:
            return False

        if p == (n - 1):
            return True

    return False


"""

performs  prime check

Input : n = number which has to be checked

Output : true or false

"""


def check_prime(n):
    if n == 2 or n == 3:
        return True

    elif n <= 1 or n % 2 == 0:
        return False

    r = n - 1

    while r % 2 == 0:
        r //= 2

    for i in range(128):
        if not miller_rabin_test(r, n):
            return False

    return True


"""

generate a prime number

"""


def give_prime(k=256):
    prime = generate_number(k)

    while not check_prime(prime):
        prime = generate_number(k)

    return prime

def give_parameters():
    return (SystemRandom().randint(0,100), SystemRandom().randint(1,100), give_prime())

"""

finds the the public key using a random starting point on the elliptic curve and a random 256 bit private key

Output : returns a Tuple of Public Key and the Original Point


"""
