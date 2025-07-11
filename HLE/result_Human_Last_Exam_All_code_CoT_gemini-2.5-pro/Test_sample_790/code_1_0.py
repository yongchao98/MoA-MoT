import math

def check(num, den, step_name):
    """Checks if a fraction is valid in the Titan architecture."""
    if num > 31 or den > 31:
        raise ValueError(
            f"Constraint violated at step '{step_name}': "
            f"Fraction {num}/{den} contains a value greater than 31."
        )
    return (num, den)

def print_fraction(f, name):
    """Prints a fraction."""
    print(f"{name} = {f[0]}/{f[1]}")

def main():
    """
    Solves the projectile problem using Titan's 5-bit fractional arithmetic.
    """
    print("--- Titan Computer Calculation Start ---")
    print("The required force F can be found using the formula: F = (2 * m * g) / cos(45°)\n")

    # Step 1: Define constants as fractions of 5-bit integers (0-31).
    print("Step 1: Define constants and approximations.")
    # Rock properties
    r = (1, 2)  # 0.5 cm
    rho = (9, 10) # 0.9 kg/cm^3
    print_fraction(r, "Radius r")
    print_fraction(rho, "Density rho")

    # Physical constants
    # For pi, 22/7 would lead to m = 33/70 (numerator > 31). We must use a less precise value.
    # Using pi=3 allows the mass calculation to stay within constraints.
    pi = (3, 1)
    print_fraction(pi, "pi (approximated)")

    # For g, a more precise value like 29/3 leads to overflow when multiplied.
    # g=10/1 is a rough approximation that keeps intermediate values valid.
    g = (10, 1)
    print_fraction(g, "g (approximated)")

    # For cos(45), we use sqrt(2) ~ 24/17 (from problem description).
    # cos(45) = sqrt(2)/2 ~ (24/17) / (2/1) = 12/17. This is a valid and accurate fraction.
    cos45 = (12, 17)
    print_fraction(cos45, "cos(45°) (approximated)")
    print("-" * 20)

    # Step 2: Calculate the projectile's mass (m).
    # m = Volume * density = (4/3 * pi * r^3) * rho
    print("Step 2: Calculate mass (m) = (4/3) * pi * r^3 * rho")
    # r^3 = (1/2)^3 = 1/8
    r_cubed = (1, 8)
    check(r_cubed[0], r_cubed[1], "r^3")
    print_fraction(r_cubed, "r^3")
    # 4/3 * pi = 4/3 * 3/1 = 4/1
    term1 = (4, 1)
    check(term1[0], term1[1], "4/3 * pi")
    print_fraction(term1, "4/3 * pi")
    # V = (4/3 * pi) * r^3 = 4/1 * 1/8 = 1/2
    V = (1, 2)
    check(V[0], V[1], "Volume")
    print_fraction(V, "Volume V")
    # m = V * rho = 1/2 * 9/10 = 9/20
    m = (9, 20)
    check(m[0], m[1], "Mass")
    print_fraction(m, "Final Mass m")
    print("-" * 20)

    # Step 3: Calculate the right-hand side of the force equation (RHS = 2*m*g).
    print("Step 3: Calculate RHS = 2 * m * g")
    # 2 * m = 2/1 * 9/20 = 9/10
    two_m = (9, 10)
    check(two_m[0], two_m[1], "2*m")
    print_fraction(two_m, "2 * m")
    # RHS = (2*m) * g = 9/10 * 10/1 = 9/1
    RHS = (9, 1)
    check(RHS[0], RHS[1], "RHS")
    print_fraction(RHS, "Final RHS")
    print("-" * 20)

    # Step 4: Solve for the Force (F).
    # F * cos(45) = RHS  =>  F = RHS / cos(45)
    print("Step 4: Solve for Force F = RHS / cos(45°)")
    print("F = (9/1) / (12/17) = (9/1) * (17/12)")
    print("This operation is INVALID because 9 * 17 = 153, which exceeds the 5-bit limit of 31.")
    print("ACTION: We must approximate the fraction 17/12 (~1.416) to proceed.")
    print("We choose the approximation 13/9 (~1.444). It is a close value and allows for simplification.")
    
    # Final calculation with approximation
    approx_inv_cos45 = (13, 9)
    # F = (9/1) * (13/9). This simplifies to 13/1.
    final_F = (13, 1)
    check(final_F[0], final_F[1], "Final Force")
    
    print("\n--- Final Titan Calculation ---")
    print("The final calculation step with the chosen approximation is:")
    print(f"F = ({RHS[0]}/{RHS[1]}) * ({approx_inv_cos45[0]}/{approx_inv_cos45[1]}) = {final_F[0]}/{final_F[1]}")
    
    # Off-line analysis to determine the error `e`
    # True F for y=10m is ~13.066 N.
    # Our calculated F is 13 N.
    # Error e = |13.0 - 13.066| = 0.066.
    # Our F=13N results in a hit at y=9.95m, which is within the target zone [9.9m, 10.1m].
    # So, the calculation is possible.
    print("\n--- Titan Computer Calculation End ---")

if __name__ == '__main__':
    main()