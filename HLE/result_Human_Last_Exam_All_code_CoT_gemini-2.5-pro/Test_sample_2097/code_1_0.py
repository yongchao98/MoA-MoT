import sympy
import math

def solve_magnetization():
    """
    Solves for the magnetization M_z(1) for n=1 and n=2,
    finds the minimum, and prints the result.
    """
    # Define symbols
    B = sympy.Symbol('B')

    # Define C(B) and its derivatives at B=1
    C = sympy.cos(sympy.pi * B / 2)
    C_prime = sympy.diff(C, B)
    C_double_prime = sympy.diff(C_prime, B)
    C_triple_prime = sympy.diff(C_double_prime, B)

    # Evaluate derivatives at B=1
    c1 = C_prime.subs(B, 1)
    c3 = C_triple_prime.subs(B, 1)

    # --- Calculation for n=1 ---
    n1 = 1
    I1 = (1 / n1**n1) * B**(4 * n1)
    I1_val = I1.subs(B, 1)
    I1_prime_val = sympy.diff(I1, B).subs(B, 1)
    I1_double_prime_val = sympy.diff(I1, B, 2).subs(B, 1)

    I0_val_n1 = (I1_prime_val - I1_val) / c1
    I0_prime_val_n1 = (I1_double_prime_val - I1_prime_val) / c1
    Mz1_n1 = I0_prime_val_n1 - I0_val_n1

    # --- Calculation for n=2 ---
    n2 = 2
    I2 = (1 / n2**n2) * B**(4 * n2)
    I2_val = I2.subs(B, 1)
    I2_prime_val = sympy.diff(I2, B).subs(B, 1)
    I2_double_prime_val = sympy.diff(I2, B, 2).subs(B, 1)
    I2_triple_prime_val = sympy.diff(I2, B, 3).subs(B, 1)
    
    # Calculate I1(1), I1'(1), I1''(1) for n=2 case
    I1_val_n2 = (I2_prime_val - I2_val) / (2 * c1)
    I1_prime_val_n2 = (I2_double_prime_val - I2_prime_val) / (2 * c1)
    
    g_val = (I2_prime_val - I2_val) / 2
    g_double_prime_val = (I2_triple_prime_val - I2_double_prime_val) / 2
    I1_double_prime_val_n2 = (g_double_prime_val * c1 - g_val * c3) / c1**2

    # Calculate I0(1), I0'(1) for n=2 case
    I0_val_n2 = (I1_prime_val_n2 - I1_val_n2) / c1
    I0_prime_val_n2 = (I1_double_prime_val_n2 - I1_prime_val_n2) / c1
    Mz1_n2 = I0_prime_val_n2 - I0_val_n2
    
    # --- Find minimum and print ---
    print(f"For n=1, M_z(1) = {Mz1_n1.evalf():.4f}")
    print(f"For n=2, M_z(1) = {Mz1_n2.evalf():.4f}")
    
    if Mz1_n1.evalf() < Mz1_n2.evalf():
        n_min = 1
        min_mz1_expr = Mz1_n1
    else:
        n_min = 2
        min_mz1_expr = Mz1_n2
        
    print(f"\nThe minimum magnetization occurs at n_min = {n_min}.")
    
    num, den = min_mz1_expr.as_numer_denom()
    
    print("\nThe final equation for the minimum magnetization is:")
    print(f"M_z(1) = {num} / {den}")
    
    print("\nThe numbers in this equation are:")
    print(f"Numerator: {num}")
    # The denominator is pi, let's print its symbolic and numeric form
    print(f"Denominator: {den} (which is approximately {den.evalf():.4f})")
    print(f"\nThe value of the minimum magnetization is {min_mz1_expr.evalf()}")


solve_magnetization()