import sys

# Titan architecture limits
MAX_VAL = 63

class TitanError(Exception):
    """Custom exception for Titan architecture errors."""
    pass

def multiply_integers(a, b):
    """Helper to multiply two integers, checking for overflow."""
    if a * b > MAX_VAL:
        raise TitanError(f"Overflow multiplying integers: {a} * {b} = {a * b}, which is > {MAX_VAL}")
    return a * b

def multiply_fractions(f1, f2):
    """Helper to multiply two fractions (tuples), checking for overflow."""
    num1, den1 = f1
    num2, den2 = f2
    if num1 * num2 > MAX_VAL or den1 * den2 > MAX_VAL:
        # This would be the point to try expansion/reduction, but for int*int it's a hard stop.
        raise TitanError(f"Overflow multiplying fractions: ({num1}/{den1}) * ({num2}/{den2})")
    return (num1 * num2, den1 * den2)

def simplify_fraction(f):
    """Simplifies a fraction by dividing numerator and denominator by their GCD."""
    num, den = f
    common = gcd(num, den)
    return (num // common, den // common)

def gcd(a, b):
    """Euclidean algorithm to find the greatest common divisor."""
    while b:
        a, b = b, a % b
    return a

def titan_calculate():
    """
    Simulates the gravitational force calculation on the Titan architecture.
    """
    print("Task: Calculate gravitational force on a 50kg probe 1km from Pandora's event horizon.")
    print("---")
    print("Step 1: Define constants as 6-bit compatible fractions.")
    
    # Using approximations that fit within 6-bit integer limits
    G_frac = (20, 3)     # Approximation of G = 6.67e-11
    rho_frac = (6, 5)    # Density = 1.2e3 kg/m^3 = 1200
    R_frac = (2, 1)      # Radius = 2e6 m
    m2_frac = (50, 1)    # Probe mass = 50 kg
    pi_frac = (3, 1)     # Using pi = 3 for simplicity to avoid early overflow
    four_thirds = (4, 3)
    
    print(f"G_frac = {G_frac[0]}/{G_frac[1]}")
    print(f"rho_frac = {rho_frac[0]}/{rho_frac[1]}")
    print(f"R_frac = {R_frac[0]}/{R_frac[1]}")
    print(f"m2_frac = {m2_frac[0]}/{m2_frac[1]}")
    print(f"pi_frac = {pi_frac[0]}/{pi_frac[1]}")
    print("---")

    print("Step 2: Justify simplifying r = Rs + d to r ~ d.")
    print("A preliminary calculation shows the Schwarzschild radius Rs is mere millimeters.")
    print("Compared to d = 1km, Rs is negligible. We simplify F = G*M*m2/r^2 to F ~ G*M*m2/d^2.")
    print("This avoids calculating Rs, but we still need to calculate M.")
    print("---")
    
    print("Step 3: Attempt to calculate the fractional part of the Force.")
    print("F_frac = G_frac * (rho_frac * 4/3 * pi_frac * R_frac^3) * m2_frac")
    print("Let's try to calculate this step-by-step, simplifying as we go to avoid overflow.")

    try:
        # Let's calculate the fractional part of the volume first: V_frac = 4/3 * pi * R^3
        print("\nCalculating Volume_frac = (4/3) * pi * R^3...")
        
        # R^3_frac
        r_squared_frac = multiply_fractions(R_frac, R_frac) # (2/1)*(2/1) = 4/1
        print(f"R_frac^2 = {r_squared_frac[0]}/{r_squared_frac[1]}")
        r_cubed_frac = multiply_fractions(r_squared_frac, R_frac) # (4/1)*(2/1) = 8/1
        print(f"R_frac^3 = {r_cubed_frac[0]}/{r_cubed_frac[1]}")

        # (4/3) * pi
        term1 = multiply_fractions(four_thirds, pi_frac) # (4/3)*(3/1) = 12/3
        print(f"(4/3) * pi = ({term1[0]}/{term1[1]})")
        term1_simple = simplify_fraction(term1) # 12/3 simplifies to 4/1
        print(f"  ... simplified to {term1_simple[0]}/{term1_simple[1]}")

        # V_frac = (4/1) * (8/1)
        v_frac = multiply_fractions(term1_simple, r_cubed_frac) # (4/1)*(8/1) = 32/1
        print(f"V_frac = (4/1) * (8/1) = {v_frac[0]}/{v_frac[1]}")
        
        # Now let's introduce other terms for the final force calculation
        # We group terms strategically to perform cancellations.
        print("\nCalculating Force_frac = (G_frac * rho_frac) * V_frac * m2_frac ...")
        
        # G_frac * rho_frac
        term2 = multiply_fractions(G_frac, rho_frac) # (20/3)*(6/5) = 120/15
        print(f"G_frac * rho_frac = (20/3) * (6/5) = {term2[0]}/{term2[1]}")
        term2_simple = simplify_fraction(term2) # 120/15 simplifies to 8/1
        print(f"  ... simplified to {term2_simple[0]}/{term2_simple[1]}")

        # Now we have to multiply the remaining terms: (8/1) * (32/1) * (50/1)
        # Let's multiply the result by V_frac
        print("\nMultiplying intermediate results: (8/1) * V_frac = (8/1) * (32/1)")
        # This is where the architecture fails. 8 * 32 = 256, which is > 63.
        # The multiply_integers helper will catch this.
        final_term1 = multiply_integers(term2_simple[0], v_frac[0])

    except TitanError as e:
        print(f"\n*** CALCULATION FAILED ***")
        print(f"Reason: {e}")
        print("\nThe Titan architecture cannot multiply 8 and 32 as their product exceeds the 6-bit limit of 63.")
        print("No further simplification or expansion tricks can resolve this integer-integer overflow.")
        print("The calculation is therefore impossible on this architecture.")
        print("\nFinal equation cannot be formed.")
        # Exit gracefully after demonstrating impossibility
        sys.exit(0)

# This part will not be reached due to the sys.exit() call upon failure.
# It's here to illustrate what the final step would have been.
# F_exp = G_exp + M_exp + m2_exp - d2_exp
# F_exp = -11 + (rho_exp + V_exp) + 0 - (d_exp*2)
# F_exp = -11 + (3 + 18) - 6 = 4
# print(f"Final Result = {final_frac[0]}/{final_frac[1]} x 10^{F_exp}")

if __name__ == '__main__':
    titan_calculate()
<<<N0>>>