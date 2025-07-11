import math

def gcd(a, b):
    """Helper function to calculate the greatest common divisor."""
    while b:
        a, b = b, a % b
    return a

def multiply_and_simplify(f1, f2):
    """
    Multiplies two fractions (n1, d1) and (n2, d2), simplifying before
    multiplication to prevent intermediate overflow.
    """
    n1, d1 = f1
    n2, d2 = f2
    
    # Simplify across the fractions
    g1 = gcd(n1, d2)
    g2 = gcd(n2, d1)
    
    num = (n1 // g1) * (n2 // g2)
    den = (d1 // g2) * (d2 // g1)
    
    return (num, den)

def solve():
    """
    Derives the calculation for the mass of the rock on the Titan architecture
    and calculates the final error.
    """
    print("Deriving the calculation for the mass of the rock.")
    print("Formula: Mass = Density * (4/3) * pi * r^3")
    print("-" * 50)

    # Step 1: Define initial values as 5-bit fractions.
    # A 5-bit integer can be from 0 to 31.
    rho = (9, 10)         # Density = 0.9
    four_thirds = (4, 3)  # Constant
    r_cubed = (1, 8)      # r = 0.5, r^3 = 0.125
    
    # For pi, we choose an approximation that is computationally feasible.
    # pi approx 22/7 leads to overflow. We use pi approx 25/8 = 3.125.
    pi_approx = (25, 8)
    
    print("Using the following 5-bit fractional approximations:")
    print(f"Density (rho)      = {rho[0]}/{rho[1]}")
    print(f"Constant           = {four_thirds[0]}/{four_thirds[1]}")
    print(f"Pi approximation   = {pi_approx[0]}/{pi_approx[1]}")
    print(f"Radius^3 (r^3)     = {r_cubed[0]}/{r_cubed[1]}")
    print("-" * 50)

    # Step 2: Perform the calculation step-by-step.
    print("Calculation steps:")
    print(f"Mass = ({rho[0]}/{rho[1]}) * ({four_thirds[0]}/{four_thirds[1]}) * ({pi_approx[0]}/{pi_approx[1]}) * ({r_cubed[0]}/{r_cubed[1]})")

    # Multiply rho and 4/3
    step1_res = multiply_and_simplify(rho, four_thirds)
    print(f"Step 1: ({rho[0]}/{rho[1]}) * ({four_thirds[0]}/{four_thirds[1]}) simplifies to {step1_res[0]}/{step1_res[1]}")
    
    # Multiply result with pi
    step2_res = multiply_and_simplify(step1_res, pi_approx)
    print(f"Step 2: ({step1_res[0]}/{step1_res[1]}) * ({pi_approx[0]}/{pi_approx[1]}) simplifies to {step2_res[0]}/{step2_res[1]}")

    # Multiply result with r^3
    n_final_unreduced = step2_res[0] * r_cubed[0]
    d_final_unreduced = step2_res[1] * r_cubed[1]
    print(f"Step 3: ({step2_res[0]}/{step2_res[1]}) * ({r_cubed[0]}/{r_cubed[1]}) results in {n_final_unreduced}/{d_final_unreduced}")

    # Step 3: Handle overflow by approximation.
    print(f"\nNOTICE: The denominator {d_final_unreduced} is > 31, which is not allowed by a 5-bit integer.")
    
    # The unrepresentable fraction is 15/32 = 0.46875.
    # We find the best representable fraction a/b where a,b <= 31.
    # The fraction 8/17 is the closest valid approximation.
    final_titan_fraction = (8, 17)
    print(f"We must approximate {n_final_unreduced}/{d_final_unreduced} with the closest valid fraction, which is {final_titan_fraction[0]}/{final_titan_fraction[1]}.")
    print("-" * 50)

    # Step 4: Final result and error calculation.
    print("Final Titan Calculation Equation:")
    print(f"Mass = ({rho[0]} / {rho[1]}) * ({four_thirds[0]} / {four_thirds[1]}) * ({pi_approx[0]} / {pi_approx[1]}) * ({r_cubed[0]} / {r_cubed[1]}) --> {step1_res[0]}/{step1_res[1]} * ({pi_approx[0]}/{pi_approx[1]}) * ({r_cubed[0]}/{r_cubed[1]}) --> {step2_res[0]}/{step2_res[1]} * ({r_cubed[0]}/{r_cubed[1]}) --> {n_final_unreduced}/{d_final_unreduced} ≈ {final_titan_fraction[0]} / {final_titan_fraction[1]}")
    
    # Calculate true value using high-precision pi
    true_mass = 0.9 * (4.0/3.0) * math.pi * (0.5**3)
    
    # Calculate Titan value from the final fraction
    titan_mass = final_titan_fraction[0] / final_titan_fraction[1]
    
    # Calculate absolute error
    abs_error = abs(true_mass - titan_mass)
    
    # Round the error to 0.001
    rounded_error = round(abs_error, 3)
    
    print(f"\nTrue mass          = {true_mass:.6f}")
    print(f"Final Titan mass   = {titan_mass:.6f}")
    print(f"Absolute error     = |{true_mass:.6f} - {titan_mass:.6f}| = {abs_error:.6f}")
    print(f"Error rounded to 0.001 = {rounded_error}")
    
    # Final answer format
    print(f"\n<<<{rounded_error}>>>")

# Execute the solver
solve()