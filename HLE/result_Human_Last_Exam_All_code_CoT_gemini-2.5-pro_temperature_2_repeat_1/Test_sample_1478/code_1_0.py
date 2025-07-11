def check_overflow(value, operation_str):
    """Checks if a value exceeds the 6-bit limit (63)."""
    if value > 63:
        print(f"COMPUTATION FAILED: The operation '{operation_str}' resulted in {value}, which exceeds the 6-bit limit of 63.")
        return True
    print(f"SUCCESS: Operation '{operation_str}' resulted in {value} (within 6-bit limit).")
    return False

def main():
    """
    Titan computer simulation to calculate the gravitational force.
    This simulation will demonstrate that the calculation is not feasible.
    """
    print("Task: Calculate gravitational force on a probe near exoplanet-turned-black-hole Pandora.")
    print("Titan Architecture Constraints: 6-bit integers (0-63) for all numerators/denominators.\n")

    # Step 1: Define constants as 6-bit friendly fractions
    # Pi is approximated as 3/1
    pi_num, pi_den = 3, 1
    # Pandora's radius R = 2000 km = 2e6 m
    R_val = 2 # Exponent is 6
    # Pandora's density rho = 1.2 t/m^3 = 1200 kg/m^3 = 12e2 kg/m^3
    rho_val = 12 # Exponent is 2
    # Volume term 4/3
    vol_num, vol_den = 4, 3

    print("--- Part 1: Calculating Mass (M) of Pandora ---")
    print(f"Formula: M = rho * (4/3) * pi * R^3")
    
    # M = (rho) * (4/3) * (pi) * (R^3)
    # Exponents are handled separately. Total exponent = 2 (from rho) + 3*6 (from R^3) = 20
    
    print("\nCalculating the fractional part of the mass...")
    # First operation: rho * (4/3) => 12 * 4/3
    term1_num = rho_val * vol_num
    if check_overflow(term1_num, f"{rho_val} * {vol_num}"):
        print("\nConclusion: Calculation cannot proceed.")
        return
        
    term1_num = term1_num // vol_den # Simplify 48 / 3 = 16
    print(f"SIMPLIFIED: ({rho_val} * {vol_num}) / {vol_den} = {term1_num}")
    
    # Second operation: result * pi => 16 * 3
    term2_num = term1_num * pi_num
    if check_overflow(term2_num, f"{term1_num} * {pi_num}"):
        print("\nConclusion: Calculation cannot proceed.")
        return

    # Third operation: result * R^3 => 48 * (2^3) = 48 * 8
    # Calculate R^3 value part first.
    R3_val = R_val ** 3
    print(f"NOTE: R^3 term value is {R_val}^3 = {R3_val} (valid).")

    final_mass_num = term2_num * R3_val
    final_op_str = f"{term2_num} * {R3_val}"
    
    print(f"\nFinal mass calculation requires: {final_op_str}")
    
    # This is the step that will fail.
    if check_overflow(final_mass_num, final_op_str):
        print("\n### FINAL CONCLUSION ###")
        print("The calculation of Pandora's mass requires multiplying 48 by 8.")
        print(f"The result is 384, which cannot be stored in a 6-bit register.")
        print("Since a key value cannot be computed without violating the hardware constraints,")
        print("the entire problem is unsolvable on the Titan architecture.")
        print("\nFinal Equation Attempted:")
        print("M_numerator = (12 * 4/3) * 3 * (2^3) = 16 * 3 * 8 = 48 * 8")


if __name__ == "__main__":
    main()
    print("\nAnswer code: N0")

<<<N0>>>