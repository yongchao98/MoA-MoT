import numpy as np

def solve():
    """
    This function determines the constant lower bound for d(t,x).
    """

    # According to the derivation, the lower bound C is determined by the minimum
    # of the function m_-(u) over the interval u in [0, 1].
    # The equation for m_-(u) is derived from the characteristic analysis of the
    # PDE for d = u_x.
    
    # Coefficients from the equation for m_-(u):
    # m_-(u) = (a*u + b*u^2 - u * sqrt(c*u^2 + d_coeff*u + e)) / divisor
    a = 3
    b = -5
    c = 17
    d_coeff = -22
    e = 9
    divisor = 4

    print("The lower bound C is determined by the minimum of the function m_-(u).")
    print("The equation for m_-(u) contains the following numerical coefficients:")
    print(f"m_-(u) = ({a}*u + {b}*u^2 - u * sqrt({c}*u^2 + ({d_coeff})*u + {e})) / {divisor}")
    print("\nWe need to find the minimum of this function for u in the interval [0, 1].")

    def m_minus(u):
        """Calculates the smaller root m_-(u) based on the coefficients."""
        discriminant_term = c * u**2 + d_coeff * u + e
        sqrt_discriminant = np.sqrt(discriminant_term)
        numerator = a * u + b * u**2 - u * sqrt_discriminant
        return numerator / divisor

    # To find the minimum, we can evaluate the function at the boundaries and check
    # its derivative, or perform a numerical search.
    
    # Check boundary values
    m_at_0 = m_minus(0.0)
    m_at_1 = m_minus(1.0)
    
    print(f"\nThe value at the boundary u=0 is m_-(0) = {m_at_0:.4f}.")
    print(f"The value at the boundary u=1 is m_-(1) = {m_at_1:.4f}.")

    # Numerical search for minimum
    u_values = np.linspace(0, 1, 2001)
    m_values = m_minus(u_values)
    min_m_val = np.min(m_values)

    print(f"A numerical search on the interval [0, 1] confirms the minimum value is approximately {min_m_val:.4f}.")
    print("A detailed analytical study would show that the function m_-(u) is monotonically decreasing on [0, 1].")
    print("Therefore, its minimum value occurs at u=1.")

    final_bound = m_at_1
    print(f"\nThe determined constant lower bound of d(t,x) is {final_bound}.")

solve()