import sympy

def solve_left_coprime_factorization():
    """
    Calculates and prints a left coprime factorization H(s) = D(s)^-1 * N(s) 
    for the transfer function:
    H(s) = [[(s-1)/(s+1), 1],
            [2/(s^2-1),   0]]
    """
    # Step 1: Define the symbolic variable and transfer function
    s = sympy.Symbol('s')
    
    # We are looking for D(s) and N(s) such that D(s)H(s) = N(s).
    # Let a generic row of D(s) be [d1(s), d2(s)].
    # The corresponding row of N(s) will be [n1(s), n2(s)].
    # The equations are:
    # n1 = d1 * (s-1)/(s+1) + d2 * 2/(s**2-1)
    # n2 = d1 * 1 + d2 * 0 = d1
    
    # Step 2: Establish conditions for polynomial solutions
    # For N(s) to be a polynomial matrix, we require n1 and n2 to be polynomials.
    # For n2 to be a polynomial, d1 must be a polynomial.
    # Let's rewrite the equation for n1 with a common denominator (s^2 - 1):
    # n1 = (d1 * (s-1)**2 + 2*d2) / (s**2 - 1)
    # For n1 to be a polynomial, the numerator must be divisible by (s^2 - 1).
    # This means the numerator must be zero at the roots of the denominator, s=1 and s=-1.
    # Let Num(s) = d1(s)*(s-1)**2 + 2*d2(s).
    # Condition at s=1: Num(1) = d1(1)*(0) + 2*d2(1) = 0  =>  d2(1) = 0.
    # Condition at s=-1: Num(-1) = d1(-1)*(-2)**2 + 2*d2(-1) = 0  =>  4*d1(-1) + 2*d2(-1) = 0  =>  2*d1(-1) + d2(-1) = 0.

    # Step 3: Find two linearly independent pairs of polynomials (d1, d2) that satisfy these conditions.
    
    # First solution pair for the first row of D(s)
    # Let's choose a simple polynomial for d11(s), e.g., d11(s) = 1.
    # The conditions on d12(s) become:
    # d12(1) = 0
    # 2*d11(-1) + d12(-1) = 2*(1) + d12(-1) = 0  =>  d12(-1) = -2.
    # A simple linear polynomial d12(s) = a*s + b can be found:
    # a + b = 0
    # -a + b = -2
    # Solving this system gives b = -1 and a = 1. So, d12(s) = s - 1.
    d11 = 1
    d12 = s - 1
    
    # Calculate the first row of N(s)
    n12 = d11
    num_n11 = d11 * (s - 1)**2 + 2 * d12
    # We simplify the expression which performs the polynomial division
    n11 = sympy.simplify(num_n11 / (s**2 - 1))

    # Second solution pair for the second row of D(s)
    # Let's choose d21(s) = s to ensure linear independence from the first row.
    # The conditions on d22(s) become:
    # d22(1) = 0
    # 2*d21(-1) + d22(-1) = 2*(-1) + d22(-1) = 0  =>  d22(-1) = 2.
    # A simple linear polynomial d22(s) = a*s + b can be found:
    # a + b = 0
    # -a + b = 2
    # Solving this system gives b = 1 and a = -1. So, d22(s) = -s + 1.
    d21 = s
    d22 = -s + 1
    
    # Calculate the second row of N(s)
    n22 = d21
    num_n21 = d21 * (s - 1)**2 + 2 * d22
    n21 = sympy.simplify(num_n21 / (s**2 - 1))

    # Step 4: Construct and print the final matrices D(s) and N(s)
    D = sympy.Matrix([
        [d11, d12],
        [d21, d22]
    ])
    N = sympy.Matrix([
        [n11, n12],
        [n21, n22]
    ])
    
    print("A left coprime factorization H(s) = D(s)^-1 * N(s) is given by:")
    print("\nD(s) =")
    sympy.pprint(D)
    print("\nN(s) =")
    sympy.pprint(N)

solve_left_coprime_factorization()