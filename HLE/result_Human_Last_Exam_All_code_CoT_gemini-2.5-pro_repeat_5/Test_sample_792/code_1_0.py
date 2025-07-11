import math

def gcd(a, b):
    """Computes the greatest common divisor of a and b."""
    return math.gcd(a, b)

def mul_frac(f1, f2, step_name=""):
    """
    Multiplies two fractions f1=(n1, d1) and f2=(n2, d2)
    under Titan constraints.
    """
    n1, d1 = f1
    n2, d2 = f2

    # Pre-simplification to keep intermediate products small
    g1 = gcd(n1, d2)
    n1_s = n1 // g1
    d2_s = d2 // g1

    g2 = gcd(n2, d1)
    n2_s = n2 // g2
    d1_s = d1 // g2

    n_res = n1_s * n2_s
    d_res = d1_s * d2_s
    
    if n_res > 31 or d_res > 31:
        # This path would fail on the Titan computer
        raise ValueError(f"Constraint violated: Operation {f1} * {f2} -> {(n_res, d_res)}")

    print(f"{step_name:<25} {n1}/{d1} * {n2}/{d2} = {n_res}/{d_res}")
    return (n_res, d_res)

def div_frac(f1, f2, step_name=""):
    """Divides fraction f1 by f2."""
    n2, d2 = f2
    return mul_frac(f1, (d2, n2), step_name)

def main():
    print("--- Titan Computation for Projectile Force ---\n")

    # Physical constants as 5-bit integer fractions
    # Using d=20m, h=10m. The ratio is what matters.
    d = (20, 1)
    h = (10, 1)
    
    # Rock properties
    # Density rho = 0.9 kg/cm^3
    rho = (9, 10) 
    # Radius r = 0.5 cm = 1/2 cm
    r = (1, 2)
    r_cubed = (r[0]**3, r[1]**3)

    # Approximations for constants
    # g ~ 9.8 m/s^2. We use 10/1 for computability.
    g_approx = (10, 1)
    # pi ~ 3.14159... We use 3/1 for computability.
    pi_approx = (3, 1)
    
    four_thirds = (4, 3)

    # --- Calculation ---
    # The order of operations is crucial to stay within constraints.
    # Formula: F = (d/h) * g * m, where m = rho * (4/3) * pi * r^3
    
    print("1. Calculating rock mass (m) in kg:")
    # Path: m = ((rho * (4/3)) * pi) * r_cubed
    try:
        m_step1 = mul_frac(rho, four_thirds, "m_step1 (rho * 4/3):")
        m_step2 = mul_frac(m_step1, pi_approx, "m_step2 ( * pi):")
        m_final = mul_frac(m_step2, r_cubed, "m_step3 ( * r^3):")
        print(f"\nCalculated mass m = {m_final[0]}/{m_final[1]} kg\n")

        print("2. Calculating force (F) in N:")
        # Path: F = ((d/h) * m) * g
        d_over_h = div_frac(d, h, "Ratio (d/h):")
        f_step1 = mul_frac(d_over_h, m_final, "F_step1 ((d/h) * m):")
        f_final = mul_frac(f_step1, g_approx, "F_step2 ( * g):")
        
        print(f"\nFinal calculated force F = {f_final[0]}/{f_final[1]} N\n")

        # --- Final Equation Output ---
        print("--- Final Equation with Titan-compatible values ---")
        # F = (d/h) * g * rho * (4/3) * pi * r^3
        print(f"Force = ({d[0]}/{d[1]} / {h[0]}/{h[1]}) * {g_approx[0]}/{g_approx[1]} * {rho[0]}/{rho[1]} * {four_thirds[0]}/{four_thirds[1]} * {pi_approx[0]}/{pi_approx[1]} * ({r[0]}/{r[1]})^3 = {f_final[0]}/{f_final[1]}")
        
        # --- Error Calculation (for verification, not part of Titan simulation) ---
        F_titan = f_final[0] / f_final[1]
        
        # High-precision calculation
        m_true = (0.9) * (4/3) * math.pi * (0.5**3)
        F_true = (20/10) * m_true * 9.8
        error = abs(F_titan - F_true)
        # The prompt asks for Y[e] in the final answer block
        # print(f"\nFor verification:")
        # print(f"High-precision force F_true ≈ {F_true:.4f} N")
        # print(f"Absolute error e = |{F_titan} - {F_true:.4f}| ≈ {error:.4f}")


    except ValueError as e:
        print(f"\nCalculation failed: {e}")
        print("\nThis combination of approximations is not computable on Titan.")

main()