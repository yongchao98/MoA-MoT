import numpy as np
from itertools import product

def get_points_in_simplex(m, D):
    """
    Generates integer points u in the D-dimensional simplex m*Q_D.
    Points satisfy u_i >= 0 and sum(u) <= m.
    """
    if D == 0:
        if m >= 0:
            yield tuple()
        return

    for u1 in range(m + 1):
        for rest in get_points_in_simplex(m - u1, D - 1):
            yield (u1,) + rest

def get_points_in_anti_simplex(m, D):
    """
    Generates integer points v in the D-dimensional simplex m*(-Q_D).
    Points satisfy v_i <= 0 and sum(v) >= -m.
    """
    for u in get_points_in_simplex(m, D):
        yield tuple(-x for x in u)

def solve_for_dimension(d):
    """
    Calculates the Ehrhart polynomial p(z) for a given dimension d,
    finds its roots, and checks the given properties.
    """
    if d < 2:
        print("Dimension d must be >= 2.")
        return

    D = d - 1  # Dimension of the cross-sections
    
    # Calculate p(n) for n = 0, 1, ..., d
    n_values = np.arange(d + 1)
    p_values = []
    
    print(f"--- Calculating for dimension d = {d} ---")
    
    for n in n_values:
        total_points = 0
        for k in range(n + 1):
            nk = n - k
            
            points_Q0 = list(get_points_in_simplex(nk, D))
            points_Q1 = list(get_points_in_anti_simplex(k, D))
            
            minkowski_sum_points = set()
            for p0 in points_Q0:
                for p1 in points_Q1:
                    minkowski_sum_points.add(tuple(p0[i] + p1[i] for i in range(D)))
            
            total_points += len(minkowski_sum_points)
        p_values.append(total_points)
        print(f"p({n}) = {total_points}")

    # Fit a polynomial of degree d to the points (n, p(n))
    coeffs = np.polyfit(n_values, p_values, d)
    
    # Round coefficients to avoid floating point inaccuracies, representing as fractions
    coeffs_frac = [np.format_float_positional(c, precision=6) for c in coeffs]

    # Print the Ehrhart polynomial
    poly_str = "p(z) = "
    for i, c in enumerate(coeffs_frac):
        power = d - i
        if abs(float(c)) > 1e-6: # Only print non-zero terms
            if i > 0 and float(c) > 0:
                poly_str += " + "
            elif float(c) < 0:
                 poly_str += " - "
                 c = str(abs(float(c)))

            poly_str += f"{c}*z^{power}" if power > 1 else (f"{c}*z" if power == 1 else f"{c}")

    print("\nEhrhart polynomial:")
    print(poly_str)

    # Find and print the roots of the polynomial
    roots = np.roots(coeffs)
    print("\nRoots of the polynomial:")
    for r in roots:
        print(f"{r:.6f}")

    # Check the answer choices
    print("\n--- Verifying Answer Choices ---")
    
    # A. Every root of p has real part -1.
    # E. Every root of p has real part -1/2.
    real_parts_A = np.allclose([r.real for r in roots], -1.0)
    real_parts_E = np.allclose([r.real for r in roots], -0.5)
    print(f"A. All roots have real part -1: {real_parts_A}")
    print(f"E. All roots have real part -1/2: {real_parts_E}")
    
    # B. Every root of p is real.
    is_real = np.allclose([r.imag for r in roots], 0)
    print(f"B. All roots are real: {is_real}")
    
    # C. The coefficients of p sum exactly d.
    # D. The coefficients of p sum exactly d choose d/2.
    sum_coeffs = p_values[1] # p(1) is the sum of coefficients
    print(f"Sum of coefficients p(1) = {sum_coeffs}")
    print(f"C. Sum is d ({d}): {sum_coeffs == d}")
    if d % 2 == 0:
        from math import comb
        d_choose_d_half = comb(d, d // 2)
        print(f"D. Sum is C(d, d/2) ({d_choose_d_half}): {sum_coeffs == d_choose_d_half}")
    else:
        print("D. d is odd, so C(d, d/2) is not an integer.")

if __name__ == '__main__':
    # We choose d=4 as a representative example beyond the trivial cases.
    solve_for_dimension(4)