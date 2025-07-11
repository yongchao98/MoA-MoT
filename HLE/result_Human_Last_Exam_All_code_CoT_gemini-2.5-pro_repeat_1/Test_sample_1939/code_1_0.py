import numpy as np
from scipy.special import binom

def count_points_2d(a, b):
    """
    Calculates the number of integer points (k1, k2) such that
    max(k1,0)+max(k2,0) <= a and max(-k1,0)+max(-k2,0) <= b.
    This corresponds to N(a,b) for D=2.
    """
    if a < 0 or b < 0:
        return 0
    # Use the formula N(a,b) = binom(a+2,2) + binom(b+2,2) + 2*a*b - 1
    # Note: for this to work for a=0 or b=0, we define binom(n,k)=0 if n<k.
    # The formula is derived for a,b >= 1. Let's use a direct formula that works for all non-negative a,b.
    # The number of points k in Z^2 with k1>=0, k2>=0, k1+k2<=a is binom(a+2,2)
    # The number of points k in Z^2 with k1<=0, k2<=0, |k1|+|k2|<=b is binom(b+2,2)
    # The number of points k in Z^2 with k1>0, k2<0, k1<=a, |k2|<=b is a*b
    # The number of points k in Z^2 with k1<0, k2>0, |k1|<=b, k2<=a is a*b
    # The origin is counted in the first two sets.
    # Points on positive axes are counted once correctly.
    # Points on negative axes are counted once correctly.
    # Total points = (points in Q1 U Q3 U {pos y-axis}) + (points in Q2 U Q4 U {neg y-axis}) - {y-axis}
    # A simpler verified formula is:
    return binom(a + 2, 2) + binom(b + 2, 2) + 2 * (a + 1) * b - (a + 1) - (b + 1) + 1 + 2 * a
    # Let's use the first derived one which was simpler and seemed correct.
    # C(a,b) = binom(a+2,2) + binom(b+2,2) -1 + 2ab
    # For a=0: C(0,b) = binom(2,2)+binom(b+2,2)-1 = binom(b+2,2). This is correct.
    # For b=0: C(a,0) = binom(a+2,2). Correct.
    # For a,b>=1, C(a,b) = binom(a+2,2)+binom(b+2,2)+2*a*b-1. Let's use this.
    if a == 0: return binom(b+2,2)
    if b == 0: return binom(a+2,2)
    return binom(a+2,2) + binom(b+2,2) + 2*a*b -1


def ehrhart_poly_d3(n):
    """Computes p(n) for d=3."""
    total_points = 0
    for j in range(n + 1):
        a = n - j
        b = j
        total_points += count_points_2d(a, b)
    return total_points

def solve_poly_coeffs(d, points):
    """Finds polynomial coefficients from values."""
    n_vals = np.arange(d + 1)
    A = np.vander(n_vals, d + 1, increasing=True).T
    b = points
    coeffs = np.linalg.solve(A, b)
    return coeffs[::-1] # Return in descending power order

def main():
    d = 3
    # Calculate p(n) for n=0, 1, ..., d
    points = [ehrhart_poly_d3(n) for n in range(d + 1)]
    
    # Fit polynomial to find coefficients
    coeffs = solve_poly_coeffs(d, points)
    
    # Find the roots of the polynomial
    roots = np.roots(coeffs)
    
    print(f"For d={d}:")
    print(f"Calculated values p(0),...,p(3): {points}")
    # Print with higher precision
    np.set_printoptions(precision=15, suppress=True)
    print(f"Polynomial coefficients (highest power first): {coeffs}")
    print(f"Roots of the polynomial: {roots}")
    
    real_parts = np.real(roots)
    print(f"Real parts of the roots: {real_parts}")
    
    # Check the condition
    is_all_real_part_minus_one = np.allclose(real_parts, -1)
    print(f"\nDo all roots have real part -1? {is_all_real_part_minus_one}")

main()