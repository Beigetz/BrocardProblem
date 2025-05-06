import sympy

# definition of symbols
p, j, m, n, k = sympy.symbols('p j m n k', integer=True)

# definition of Kronecker Delta function
def delta(i, q):
    return sympy.KroneckerDelta(i, q)

# definition of m-th multifactorial(1 if n < 1)
def multifactorial(val, step):
    if step <= 0:
        pass
    if val < 1:
        return 1
    result = 1
    current = val
    while current >= 1:
        result *= current
        current -= step
    return result

# definition of b(p, j, m, n)
def b_func(p, j, m, n):
    product_term = \
    sympy.product(((2**(j-1) + 2*k + 1)*m + n) *\
    ((2*k + 2)*m + n), (k, 0, 2**(j-2) - 1))
    return product_term**p

# definition of Delta(p, j, m, n)
def delta_func(p, j, m, n):
    term1 = \
    (sympy.factorial(2**j * m + n) / sympy.factorial(m)**\
    (sympy.floor((2**j * m + n)/m)))**p / b_func(p, j, m, n)
    term2 = delta(m, 1) * delta(n, -1)
    term3 = multifactorial(n, m)**p if n >= 1 else 1
    if n == -1:
        return (multifactorial(2**j * m - 1, m)**p /\
        b_func(p, j, m, n)) * delta(m, 1)
    else:
        return term1 * term2 * term3
# definition of a(p, j, m, n)
def a_func(p, j, m, n):
    product_term = \
    sympy.product(((2**(j-1) + 2*k + 2)*m + n) *\
    ((2*k + 1)*m + n), (k, 0, 2**(j-2) - 1))
    n_fact_p = multifactorial(n, m)**p if n >= 1 else 1
    return product_term**p * n_fact_p +\
    delta_func(p, j, m, n)

# definition of c(p, j, m, n)
def c_func(p, j, m, n):
    return (a_func(p, j, m, n) - b_func(p, j, m, n)) / 2

# definition of square root of the right-hand side of formula
def kf(p,j,m,n):
    return a_func(p,j,m,n)-c_func(p,j,m,n)

# definition of Absolute value of c(p,j,m,n)
def cf(p,j,m,n):
    return sympy.Abs(c_func(p,j,m,n))

# definition of output of the formula under conditions
# p-th power of m-th multifuctorial of
# Brocard-Ramanujan-type Diophantine equation
def breq(p,j,m,n):
    # Check integer or not
    if (isinstance(p, int) or \
    ( p.is_integer() and isinstance(p,float)))\
    and (isinstance(j, int) or \
    ( j.is_integer() and isinstance(j,float)))\
    and (isinstance(m, int) or \
    ( m.is_integer() and isinstance(m,float)))\
    and (isinstance(n, int) or \
    ( n.is_integer() and isinstance(n,float))):

        # if each value is integer, converting float type to integer
        p=int(p); j=int(j); m=int(m); n=int(n)
        if p<=0 or j<=1 or m<=0 or n<=-2:
            print("Enter an integer that satisfies the condition.")
            pass

        # To ensure integer solutions to the Diophantine
        # equation, we exclude
        # the case (p,j,m,n)=(p,2,1,-1), i.e. (3!)^p.
        if j==2 and m==1 and n==-1:
            pass
        # Results are output only when the condition
        # for the equality to hold, i.e., the left-hand
        # side minus the right-hand side equals zero,
        # is satisfied.
        elif multifactorial(((2**j)*m+n),m)**p+\
        (cf(p,j,m,n))**2-kf(p,j,m,n)**2==0:
            if p==1:
                if m==1:
                    if cf(p,j,m,n)==1:
                        print(f"p:{p}, j:{j}, m:{m}, n:{n}")
                        print("solution: The Brocard number")
                        print(f"({(2**j*m+n)}, {kf(p,j,m,n)})")
                        print(f"{2**j*m+n}!+{cf(p,j,m,n)}"\
                        +f"={kf(p,j,m,n)}^2")
                    else:
                        pass
                else:
                    print(f"p:{p}, j:{j}, m:{m}, n:{n}")
                    print("solution:")
                    print(f"{2**j*m+n}!_{m}+({cf(p,j,m,n)})^2")
                    print(f"=({kf(p,j,m,n)})^2")
            elif p != 1:
                if m==1:
                    print(f"p:{p}, j:{j}, m:{m}, n:{n}")
                    print("solution:")
                    print(f"({(2**j*m+n)}!)^{p}+({cf(p,j,m,n)})^2")
                    print(f"=({kf(p,j,m,n)})^2")
                else:
                    print(f"p:{p}, j:{j}, m:{m}, n:{n}")
                    print("solution:")
                    print(f"(({(2**j*m+n)})!_{m})^{p}"\
                    +f"({cf(p,j,m,n)})^2")
                    print(f"=({kf(p,j,m,n)})^2")
        else:
            pass
    else:
        print("Enter an integer that satisfies the condition.")

    # This code is written with an emphasis on clarity
    # and is not optimized for high-speed implementation.
    # Checking with loops:
    # p >= 1. this case 1 <= p <= 3
    for p in range(1,4):
        # j >= 2.  this case 2 <= j <= 3
        for j in range(2,4):
            # m >= 1.  this case 1 <= m <= 3
            for m in range(1,4):
                # n >= -1.  this case -1 <= n <= 9
                for n in range(-1,10):
                    breq(p,j,m,n)
