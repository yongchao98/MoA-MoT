import math

def check_constraint(n, d):
    """Checks if a fraction respects the 5-bit constraint."""
    if n > 31 or d > 31:
        raise ValueError(f"Constraint violated: fraction {n}/{d}. Components must be <= 31.")
    return True

def simplify_fraction(n, d):
    """Simplifies a fraction using the greatest common divisor."""
    if d == 0:
        raise ValueError("Denominator cannot be zero.")
    common = math.gcd(n, d)
    return n // common, d // common

def multiply_fractions(f1, f2):
    """
    Multiplies two fractions with simplification to stay within Titan constraints.
    This special multiplication function simplifies before multiplying the numerators
    and denominators to avoid large intermediate numbers that would violate the constraints.
    """
    n1, d1 = f1
    n2, d2 = f2
    
    # Cross-simplify before multiplication
    g1 = math.gcd(n1, d2)
    g2 = math.gcd(n2, d1)
    
    new_n = (n1 // g1) * (n2 // g2)
    new_d = (d1 // g2) * (d2 // g1)
    
    # The result of multiplication itself must not violate the constraint
    check_constraint(new_n, new_d)
    
    return new_n, new_d

def divide_fractions(f1, f2):
    """Divides two fractions by multiplying by the reciprocal."""
    n2, d2 = f2
    if n2 == 0:
        raise ValueError("Cannot divide by zero fraction.")
    return multiply_fractions(f1, (d2, n2))

def power_fraction(f, p):
    """Calculates the power of a fraction."""
    n, d = f
    res_n, res_d = n**p, d**p
    check_constraint(res_n, res_d)
    return res_n, res_d

# --- Main Calculation ---
try:
    # 1. Define constants and approximations as fractions (numerator, denominator)
    # Problem values
    radius = (1, 2)       # 0.5 cm
    density = (9, 10)     # 0.9 kg/cm^3
    height = (10, 1)      # 10 m
    distance = (20, 1)    # 20 m
    
    # Universal constants approximated as fractions
    four_thirds = (4, 3)
    pi_approx = (3, 1)    # approx for pi ≈ 3.14159
    g_approx = (10, 1)    # approx for g ≈ 9.8 m/s^2

    # The equation for the force F is: F = (d/h) * density * (4/3) * pi * r^3 * g
    # We will compute this step-by-step to respect Titan's computational rules.
    print("Deriving the force F for the Titan computer.")
    print("Formula: F = (distance / height) * density * (4/3) * (radius^3) * pi * g")
    print(f"Using values: distance={distance[0]}/{distance[1]}, height={height[0]}/{height[1]}, density={density[0]}/{density[1]}, pi={pi_approx[0]}/{pi_approx[1]}, g={g_approx[0]}/{g_approx[1]}, radius={radius[0]}/{radius[1]}")
    print("-" * 30)

    # Perform the calculation in an order that respects the constraints
    print("Calculation steps:")
    # Step 1: radius^3
    radius_cubed = power_fraction(radius, 3)
    print(f"radius^3 = ({radius[0]}/{radius[1]})^3 = {radius_cubed[0]}/{radius_cubed[1]}")

    # The rest of the coefficients: C = (d/h) * density * (4/3) * g * pi
    # Grouping them in a specific order is crucial to stay within constraints.
    # C1 = d/h
    term1 = divide_fractions(distance, height)
    print(f"(distance/height) = ({distance[0]}/{distance[1]}) / ({height[0]}/{height[1]}) = {term1[0]}/{term1[1]}")

    # C2 = C1 * g
    term2 = multiply_fractions(term1, g_approx)
    print(f"-> multiplied by g ({g_approx[0]}/{g_approx[1]}) = {term2[0]}/{term2[1]}")
    
    # C3 = C2 * density
    term3 = multiply_fractions(term2, density)
    print(f"-> multiplied by density ({density[0]}/{density[1]}) = {term3[0]}/{term3[1]}")

    # C4 = C3 * 4/3
    term4 = multiply_fractions(term3, four_thirds)
    print(f"-> multiplied by (4/3) = {term4[0]}/{term4[1]}")
    
    # C5 = C4 * pi
    term5 = multiply_fractions(term4, pi_approx)
    print(f"-> multiplied by pi ({pi_approx[0]}/{pi_approx[1]}) = {term5[0]}/{term5[1]}")
    
    # Final step: F = C5 * r^3
    final_F_frac = multiply_fractions(term5, radius_cubed)
    print(f"-> multiplied by radius^3 ({radius_cubed[0]}/{radius_cubed[1]}) = {final_F_frac[0]}/{final_F_frac[1]}")
    
    print("-" * 30)
    
    # Final calculated force
    force_titan = final_F_frac[0] / final_F_frac[1]
    print(f"Final calculated Force (F) by Titan: {final_F_frac[0]}/{final_F_frac[1]} = {force_titan} N")
    
    # For comparison, calculate the "true" force for d=20m using float values
    r_true, d_true, h_true, dens_true, g_true, pi_true = 0.005, 20.0, 10.0, 900000, 9.8, math.pi
    # m = density * volume, need to convert units properly
    # r=0.5cm=0.005m, density=0.9kg/cm3=900,000 kg/m3
    m_true = dens_true * (4/3) * pi_true * (r_true**3)
    force_true = (d_true / h_true) * m_true * g_true
    print(f"True required force for d=20m: {force_true:.4f} N")

    # Check if the calculated force hits the target
    d_actual = (h_true * force_titan) / (m_true * g_true)
    print(f"The actual landing distance with Titan's calculated force would be: {d_actual:.4f} m")

    if 19 <= d_actual <= 21:
        print("This is within the 19-21m range, so the calculation is a success.")
        error = abs(force_titan - force_true)
        print(f"Absolute error in calculated force: {error:.4f}")
        # Final answer format
        final_answer = f"Y[{error:.4f}]"
    else:
        print("This is outside the 19-21m range, so the calculation failed.")
        final_answer = "N0"

except ValueError as e:
    print(f"An error occurred: {e}")
    final_answer = "N0"
    
print("<<<" + final_answer + ">>>")
