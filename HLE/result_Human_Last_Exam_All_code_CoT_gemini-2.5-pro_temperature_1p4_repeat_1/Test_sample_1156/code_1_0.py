import numpy as np
from scipy.special import zeta

def solve_density_approximation():
    """
    Calculates the parameters for an approximate invariant density rho(x) = C * (x**(-1/2) + b).
    The parameters are found by enforcing the Perron-Frobenius equation at y=1.
    """
    
    # Zeta function values
    # zeta(2) is pi^2/6
    zeta_2 = np.pi**2 / 6
    zeta_3 = zeta(3)
    
    # Calculate b from the condition at y=1
    b_numerator = 2 * zeta_2 - 3
    b_denominator = 3 - 2 * zeta_3
    b = b_numerator / b_denominator
    
    # Calculate the normalization constant C
    C = 1 / (2 + b)
    
    # The coefficients of the final equation rho(x) = C * (x**p1 + b)
    p1 = -0.5
    
    print("The approximate normalised density of the invariant measure is rho(x) = C * (x**p + b)")
    print("where the parameters are:")
    print(f"C = {C}")
    print(f"p = {p1}")
    print(f"b = {b}")
    print("\nSo the equation is:")
    print(f"rho(x) = {C} * (x**({p1}) + ({b}))")

solve_density_approximation()