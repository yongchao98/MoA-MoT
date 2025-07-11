import numpy as np

def solve():
    """
    This script calculates the infimum and supremum of the expression 
    |P_n(xi)| * (a_n^2 + b_n^2 + c_n^2) based on a theoretical formula.
    """
    
    # Define the polynomial f(x) = x^3 + x^2 + x - 1
    f_coeffs = [1, 1, 1, -1]
    
    # Find the roots of the polynomial
    roots = np.roots(f_coeffs)
    
    # The roots are xi (real), xi_1 (complex), and xi_2 (complex conjugate)
    # The real root is approximately 0.54
    real_root_mask = np.isreal(roots)
    xi = roots[real_root_mask][0].real
    
    # The complex roots are the other two
    complex_roots = roots[~real_root_mask]
    xi_1 = complex_roots[0]
    
    # Discriminant of the polynomial x^3+px^2+qx+r=0 is p^2q^2-4q^3-4p^3r-27r^2+18pqr
    # For f(x) = x^3+x^2+x-1, p=1, q=1, r=-1
    p, q, r = 1, 1, -1
    discriminant = p**2 * q**2 - 4 * q**3 - 4 * p**3 * r - 27 * r**2 + 18 * p * q * r
    # So, Delta = -44
    
    # Define the derivative f'(x) = 3x^2 + 2x + 1
    df = lambda x: 3*x**2 + 2*x + 1
    
    # Evaluate the derivative at the roots
    df_xi = df(xi)
    df_xi_1 = df(xi_1)
    
    # Get the magnitudes
    abs_df_xi = np.abs(df_xi)
    abs_df_xi_1 = np.abs(df_xi_1)
    
    # The theoretical formula for the endpoints of the limit interval
    # is (|Delta| / (2 * |f'(xi_1)|^2)) * [|f'(xi_1)|/|f'(xi)| +/- |f'(xi_1)/f'(xi) - 1|]
    
    abs_delta = np.abs(discriminant)
    
    term1 = abs_delta / (2 * abs_df_xi_1**2)
    
    ratio_mag = abs_df_xi_1 / abs_df_xi
    ratio_complex_mag = np.abs(df_xi_1 / df_xi - 1)
    
    inf_val = term1 * (ratio_mag - ratio_complex_mag)
    sup_val = term1 * (ratio_mag + ratio_complex_mag)
    
    print(f"The discriminant Delta is: {discriminant}")
    print(f"The real root xi is approximately: {xi:.6f}")
    print(f"The complex root xi_1 is approximately: {xi_1.real:.6f} + {xi_1.imag:.6f}i")
    print(f"The derivative f'(x) = 3x^2 + 2x + 1")
    print(f"|f'(xi)| is approximately: {abs_df_xi:.6f}")
    print(f"|f'(xi_1)| is approximately: {abs_df_xi_1:.6f}")
    print("-" * 30)
    print(f"The calculated infimum is: {inf_val}")
    print(f"The calculated supremum is: {sup_val}")

solve()