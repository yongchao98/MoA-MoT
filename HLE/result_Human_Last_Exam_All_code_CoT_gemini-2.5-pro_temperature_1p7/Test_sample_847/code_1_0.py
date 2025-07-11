import numpy as np

def solve_minimal_cost():
    """
    Calculates the minimal number of coins per bit of information needed to sort an array.
    This is the core of solving the problem for large n.
    """
    
    # The optimal probability 'p' for a "yes" answer minimizes the cost per bit.
    # It is the real root of the cubic equation p^3 + p - 1 = 0.
    # We define the coefficients of the polynomial: 1*p^3 + 0*p^2 + 1*p - 1.
    coeffs = [1, 0, 1, -1]
    
    # Find the roots of the polynomial.
    roots = np.roots(coeffs)
    
    # The equation has one real root and two complex conjugate roots. We need the real one.
    p_optimal = roots[np.isreal(roots)].real[0]
    
    # The minimal cost per bit of information is given by the formula -log_p(2),
    # which can be calculated as -ln(2) / ln(p).
    # where p is the optimal probability we just found.
    cost_per_bit = -np.log(2) / np.log(p_optimal)
    
    print("The optimal probability 'p' is the real root of the equation p^3 + p - 1 = 0.")
    print(f"p = {p_optimal:.7f}")
    print("\nThe minimal cost per bit is calculated using the formula: -log_p(2)")
    print(f"Minimal cost = -log_({p_optimal:.7f})(2) = {cost_per_bit:.7f}")
    
    print("\nRounded to 3 decimal places, the minimal number of coins (cost per bit) is:")
    print(f"{cost_per_bit:.3f}")

solve_minimal_cost()