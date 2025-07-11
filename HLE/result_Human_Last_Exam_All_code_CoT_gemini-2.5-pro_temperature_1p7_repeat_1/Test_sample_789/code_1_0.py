import math

def gcd(a, b):
    """Computes the greatest common divisor of two integers."""
    while b:
        a, b = b, a % b
    return a

def find_best_approximation(target_num, target_den):
    """
    Finds the best fraction a/b with a, b <= 31 to approximate target_num / target_den.
    """
    best_fraction = (0, 1)
    min_error = float('inf')
    
    target_value = target_num / target_den
    
    # Iterate through all possible denominators (1 to 31)
    for b in range(1, 32):
        # Find the closest integer numerator for this denominator
        ideal_a = round(b * target_value)
        # Ensure numerator is also within the 5-bit limit
        if 0 <= ideal_a <= 31:
            a = int(ideal_a)
            if a == 0 and b == 1: continue # Skip 0/1 unless target is 0

            current_error = abs(target_value - a / b)
            if current_error < min_error:
                min_error = current_error
                best_fraction = (a, b)
    
    return best_fraction

def main():
    """
    Performs the calculation to find the mass of the rock on the Titan architecture.
    """
    print("Deriving the calculation for the mass of the rock.")
    print("The formula is: Mass = Density * Volume, where Volume = (4/3) * pi * r^3")
    print("-" * 50)

    # Step 1: Define initial values as fractions with 5-bit integers
    density_n, density_d = 9, 10
    radius_n, radius_d = 1, 2
    four_thirds_n, four_thirds_d = 4, 3
    pi_n, pi_d = 22, 7  # Using the approximation pi ~= 22/7
    
    print(f"Initial values as 5-bit fractions:")
    print(f"Density (rho) = {density_n}/{density_d}")
    print(f"Radius (r) = {radius_n}/{radius_d}")
    print(f"Pi approximation = {pi_n}/{pi_d}")
    print("-" * 50)

    # Step 2: Show the full expression
    print("The full calculation is:")
    print(f"Mass = ({density_n}/{density_d}) * ({four_thirds_n}/{four_thirds_d}) * ({pi_n}/{pi_d}) * ({radius_n}/{radius_d})^3")
    print("-" * 50)

    # Step 3: Perform calculation step-by-step
    print("Step-by-step calculation adhering to 5-bit constraints:")

    # Calculate r^3
    r_cubed_n, r_cubed_d = radius_n**3, radius_d**3
    print(f"1. Calculate r^3: ({radius_n}/{radius_d})^3 = {r_cubed_n}/{r_ cubed_d}")

    # Combine pi * r^3
    # M = (9/10) * (4/3) * (22/7) * (1/8)
    # To avoid overflow, we reorder: M = (9/10) * [ ( (22/7) * (1/8) ) * (4/3) ]
    print("2. Combine terms carefully to avoid intermediate overflow.")
    
    # op1 = pi * r^3 = 22/7 * 1/8
    common_divisor_1 = gcd(pi_n, r_cubed_d)
    op1_n = (pi_n // common_divisor_1) * r_cubed_n
    op1_d = pi_d * (r_cubed_d // common_divisor_1)
    print(f"   a. Pi * r^3 = {pi_n}/{pi_d} * {r_cubed_n}/{r_cubed_d} -> simplify by 2 -> {op1_n}/{op1_d}")
    
    # op2 = (4/3) * op1 = 4/3 * 11/28
    common_divisor_2 = gcd(four_thirds_n, op1_d)
    op2_n = (four_thirds_n // common_divisor_2) * op1_n
    op2_d = four_thirds_d * (op1_d // common_divisor_2)
    print(f"   b. (4/3) * (result from 2a) = {four_thirds_n}/{four_thirds_d} * {op1_n}/{op1_d} -> simplify by 4 -> {op2_n}/{op2_d}")
    
    # op3 = density * op2 = 9/10 * 11/21
    print(f"   c. Density * (result from 2b) = {density_n}/{density_d} * {op2_n}/{op2_d}")
    
    common_divisor_num = gcd(density_n, op2_d)
    common_divisor_den = gcd(density_d, op2_n) # Not possible in this case
    
    final_ideal_n = (density_n // common_divisor_num) * (op2_n // common_divisor_den)
    final_ideal_d = (density_d // common_divisor_den) * (op2_d // common_divisor_num)
    
    print(f"      This simplifies to {final_ideal_n}/{final_ideal_d}.")
    print(f"      Problem: The numerator ({final_ideal_n}) and denominator ({final_ideal_d}) are larger than 31 and cannot be stored.")
    print("-" * 50)
    
    # Step 4: Find the best 5-bit approximation for the unrepresentable result
    print(f"Step 4: Find the best 5-bit fraction to approximate the result {final_ideal_n}/{final_ideal_d}.")
    
    calc_n, calc_d = find_best_approximation(final_ideal_n, final_ideal_d)
    print(f"The best representable fraction is {calc_n}/{calc_d}.")
    print("-" * 50)

    # Step 5: Final results
    print("Final Calculation Summary:")
    print(f"Mass = ({density_n}/{density_d}) * ({four_thirds_n}/{four_thirds_d}) * ({pi_n}/{pi_d}) * ({radius_n}/{radius_d})^3  (approximated to) ==> {calc_n}/{calc_d}")
    
    # Calculate true mass and absolute error
    true_mass = (9/10) * (4/3) * math.pi * (1/2)**3
    calculated_mass_val = calc_n / calc_d
    absolute_error = abs(true_mass - calculated_mass_val)
    
    print(f"\nTrue mass = {true_mass:.6f}")
    print(f"Calculated mass on Titan = {calculated_mass_val:.6f}")
    print(f"Smallest absolute error = {absolute_error:.6f}")
    
    rounded_error = round(absolute_error, 3)
    print(f"Smallest absolute error rounded to 0.001 = {rounded_error}")
    
    # Final answer in requested format
    print(f"\n<<<{rounded_error}>>>")


if __name__ == "__main__":
    main()