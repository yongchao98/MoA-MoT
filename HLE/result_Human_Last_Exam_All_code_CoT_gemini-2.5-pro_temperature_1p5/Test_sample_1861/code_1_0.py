import sympy

def solve():
    """
    This script verifies the counterexample for the case M = R^2.
    It shows that the 1-form eta = exp(x)*dy has a non-zero exterior derivative
    and is invariant under a transitive group of diffeomorphisms.
    This demonstrates that d(eta) is not necessarily zero for R^2.
    """

    # Define symbols for coordinates and transformation parameters
    x, y, c, d = sympy.symbols('x y c d')

    # Define the 1-form eta. In this representation, it's a list [f, g] for f*dx + g*dy
    # Our eta = e^x*dy, so the dx component is 0 and dy component is exp(x)
    eta_dx = 0
    eta_dy = sympy.exp(x)
    print(f"Let eta = ({eta_dx}) dx + ({eta_dy}) dy")

    # Compute the exterior derivative d(eta) = (dg/dx - df/dy) dx^dy
    # where eta = f dx + g dy
    d_eta_coeff = sympy.diff(eta_dy, x) - sympy.diff(eta_dx, y)
    print(f"The exterior derivative d(eta) is ({d_eta_coeff}) dx^dy.")
    print(f"Since {d_eta_coeff} = exp(x) is never zero, d(eta) is not identically zero.")
    print("-" * 20)

    # Define a transitive transformation F(x, y) = (x', y')
    # Let F(x,y) = (x - ln(c), c*y + d). This transformation is transitive on R^2 for c > 0.
    x_prime = x - sympy.log(c)
    y_prime = c * y + d
    print(f"Consider the transitive diffeomorphism F(x,y) = (x', y') where:")
    print(f"x' = {x_prime}")
    print(f"y' = {y_prime}")
    print("-" * 20)


    # To compute the pullback F*(eta), we use the formula:
    # F*(f dx + g dy) = (f(F) * dx'/dx + g(F) * dy'/dx) dx + (f(F) * dx'/dy + g(F) * dy'/dy) dy
    
    # Substitute x' and y' into eta's coefficients
    eta_dx_F = eta_dx.subs([(x, x_prime), (y, y_prime)])
    eta_dy_F = eta_dy.subs([(x, x_prime), (y, y_prime)])

    # Compute Jacobian matrix elements of F
    dx_prime_dx = sympy.diff(x_prime, x)
    dx_prime_dy = sympy.diff(x_prime, y)
    dy_prime_dx = sympy.diff(y_prime, x)
    dy_prime_dy = sympy.diff(y_prime, y)
    
    # Compute the coefficients of the pullback form F*(eta)
    pullback_eta_dx = eta_dx_F * dx_prime_dx + eta_dy_F * dy_prime_dx
    pullback_eta_dy = eta_dx_F * dx_prime_dy + eta_dy_F * dy_prime_dy

    # Simplify the results
    pullback_eta_dx_simple = sympy.simplify(pullback_eta_dx)
    pullback_eta_dy_simple = sympy.simplify(pullback_eta_dy)
    
    print("Now we compute the pullback F*(eta):")
    print(f"The dx component of F*(eta) is: {pullback_eta_dx_simple}")
    print(f"The dy component of F*(eta) is: {pullback_eta_dy_simple}")

    # Check if F*(eta) == eta
    if pullback_eta_dx_simple == eta_dx and pullback_eta_dy_simple == eta_dy:
        print("\nConclusion: F*(eta) = eta. The form is invariant under the transformation.")
        print("Since we found a valid 1-form eta on R^2 with d(eta) != 0 that satisfies the condition,")
        print("it is NOT necessary for d(eta) to be 0 on R^2.")
    else:
        print("\nConclusion: F*(eta) != eta. The calculation is incorrect or the form is not invariant.")
        
solve()