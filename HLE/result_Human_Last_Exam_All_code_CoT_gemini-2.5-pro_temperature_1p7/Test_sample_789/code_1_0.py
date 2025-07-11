import math

def calculate_mass_and_error():
    """
    This function performs the entire calculation as per the Titan 5-bit computer rules
    and then calculates the final error.
    """

    # Step 1 & 2: Define constants as 5-bit integer fractions
    # Density (rho) = 0.9 kg/cm^3 = 9/10
    rho_num, rho_den = 9, 10
    
    # Radius (r) = 0.5 cm = 1/2
    r_num, r_den = 1, 2
    
    # Constant 4/3
    four_thirds_num, four_thirds_den = 4, 3

    print("Titan Calculation for Mass of the Rock")
    print("========================================")
    print(f"Formula: Mass = Density * (4/3) * pi * r^3")
    print(f"Inputs:")
    print(f"  Density (rho) = {rho_num}/{rho_den}")
    print(f"  Radius (r) = {r_num}/{r_den}")
    print(f"  Constant = {four_thirds_num}/{four_thirds_den}")

    # Step 3: Choose an approximation for pi
    # The most accurate simple fraction, 22/7, leads to M = 33/70, which overflows the 5-bit registers.
    # We choose a less precise but computable approximation, pi = 28/9, to allow for simplification.
    pi_num, pi_den = 28, 9
    print(f"  Chosen pi = {pi_num}/{pi_den}")
    
    # Calculate r^3
    r3_num = r_num**3
    r3_den = r_den**3
    print(f"Step 1: Calculate r^3 = ({r_num}/{r_den})^3 = {r3_num}/{r3_den}")

    # Step 4: Full calculation by simplifying a pool of numerators and denominators
    print("Step 2: Combine all fractions into a single expression.")
    # M = (9/10) * (4/3) * (28/9) * (1/8)
    numerators = [rho_num, four_thirds_num, pi_num, r3_num]
    denominators = [rho_den, four_thirds_den, pi_den, r3_den]
    print(f"  Mass = ({numerators[0]}*{numerators[1]}*{numerators[2]}*{numerators[3]}) / ({denominators[0]}*{denominators[1]}*{denominators[2]}*{denominators[3]})")

    # To stay within 5-bit limits, we simplify before multiplying the pools.
    # We can think of it as a pool of numerators {9, 4, 28, 1} and denominators {10, 3, 9, 8}.
    print("Step 3: Simplify the expression by cancelling common factors.")
    print("  Original Numerators: {9, 4, 28, 1}, Denominators: {10, 3, 9, 8}")
    
    # Cancel 9 from num and den
    print("  - Cancel factor 9 from numerator and denominator.")
    # New pools: Num: {4, 28, 1}, Den: {10, 3, 8}
    
    # Cancel 4 from num and 8 from den (leaves 2 in den)
    print("  - Cancel factor 4 from numerator and denominator (4 in num, 8 in den leaves 2).")
    # New pools: Num: {28, 1}, Den: {10, 3, 2}

    # Cancel 28 from num and 2 from den (leaves 14 in num)
    print("  - Cancel factor 2 from numerator and denominator (28 in num, 2 in den leaves 14).")
    # New pools: Num: {14, 1}, Den: {10, 3}

    final_num = 14
    final_den = 30
    print(f"Step 4: Multiply the remaining factors: {final_num}/{final_den}")

    # Simplify the final fraction
    common_divisor = math.gcd(final_num, final_den)
    result_num = final_num // common_divisor
    result_den = final_den // common_divisor
    print(f"Step 5: Reduce the final fraction: {result_num}/{result_den}")
    
    # The result of the calculation
    calculated_mass = result_num / result_den
    print("\nFinal Titan Calculation Result:")
    print(f"Mass = {result_num} / {result_den} kg (approx. {calculated_mass:.4f} kg)")

    # Step 5: Calculate the error
    print("\nError Calculation")
    # True value calculation using high-precision pi
    true_mass = (9/10) * (4/3) * math.pi * (1/2)**3
    print(f"  True mass (using high-precision pi) is approx. {true_mass:.7f} kg")
    
    absolute_error = abs(true_mass - calculated_mass)
    print(f"  Absolute error = |{true_mass:.7f} - {calculated_mass:.7f}| = {absolute_error:.7f}")
    
    # Step 6: Round error to 0.001
    rounded_error = round(absolute_error, 3)
    print(f"  The error rounded to 0.001 is: {rounded_error}")

    return rounded_error

# Execute the calculation and get the final answer
final_error = calculate_mass_and_error()

# Print the final answer in the required format
print(f"\n<<<{final_error}>>>")
