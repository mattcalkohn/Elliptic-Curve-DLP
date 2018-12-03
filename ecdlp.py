###############################################
#                                             #
#  Elliptic Curve Discrete Logarithm Problem  #
#                  Matt Kohn                  #
#                                             #
###############################################

# Solves the DLP for elliptic curves of the form:
#       Y**2 = X**3 + A*X + B
# in finite fields of prime order (F = {0, 1, ..., p-1})
#  for some prime p. Namely, solves the problem of
#  finding an integer n such that nP = Q for points
#  P, Q on a given elliptic curve, if such an n exists
#  The point O (the point at infinity) is defined as (None, None)

# computes the y**2 value of the elliptic curve
def f(x, a, b, p):
    return (x**3 + a*x + b) % p

# computes the gcd of a and b
def gcd(a, b):
    if a == 0:
        return (b, 0, 1)
    else:
        g, y, x = gcd(b % a, a)
        return (g, x - (b // a) * y, y)

# checks for the existence of a modular inverse to a modulo p
def modular_inverse(a, p):
    g, x, y = gcd(a, p)
    if g != 1:
        return None
    else:
        return x % p

# checks to see if there exists an x to den * x = num (mod p)
def modular_division(num, den, p):
    x = 0

    m = modular_inverse(den, p)
    if m != None:
        x = (m * num) % p

    return x

# computes the slope of the line secant to p and q
def slope(p, a, b, px, py, qx, qy):
    # variables for numerator and denominator of slope
    num = 0
    den = 0

    if px == qx and py == qy:
        # tangent line
        num = (3 * (px ** 2) + a) % p 
        den = (2 * py) % p
    elif px == qx:
        # vertical line
        return None
    else:
        # secant line
        num = (py - qy) % p
        den = (px - qx) % p

    if num % den != 0:
        # check if can do  modular inverse
        return modular_division(num, den, p)
    else:
        return (num / den) % p

# adds points p and q for the elliptic curve defined by a, b, p
def add(p, a, b, px, py, qx, qy):
    # check if either point is O
    if px == None:
        return (qx, qy)
    elif qx == None:
        return (px, py)

    dydx = slope(p, a, b, px, py, qx, qy)
    # check if p + q = O
    if dydx == None:
        return None, None

    x = (dydx ** 2 - px - qx) % p
    y = (dydx * (px - x) - py) % p

    return (x, y)

# computes nP for some integer n and a point P
def multiply(p, a, b, n, px, py):
    if n == 0:
        return None, None

    n0 = 1
    tempx = px
    tempy = py
    while n0 < n:
        tempx, tempy = add(p, a, b, tempx, tempy, px, py)
        n0 += 1

    return (tempx, tempy)

# naively solves the discrete log problem for the elliptic
#  curve defined by y**2 = x**3 + a*x + b by finding some n 
#  such that nP = Q where P = (px, py) and Q = (qx, qy)
def ecdlp(p, a, b, px, py, qx, qy):
    # check if P = O
    if px == None:
        return -1

    # check to make sure both P and Q are on the curve
    if f(px, a, b, p) != (py**2 % p) or f(qx, a, b, p) != (qy**2 % p):
        return -1

    n = 1
    tempx = px
    tempy = py
    while (tempx != qx or tempy != qy) and n < p:
        tempx, tempy = add(p, a, b, tempx, tempy, px, py)
        n += 1

    return n
