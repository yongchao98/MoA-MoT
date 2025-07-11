import math

def find_best_approximation(target_value, limit):
    """
    Finds the best fractional approximation n/d for a target value,
    where n and d are less than or equal to the given limit.
    """
    best_fraction = (0, 1)
    min_error = float('inf')

    for d in range(1, limit + 1):
        # Find the best numerator for the current denominator
        # n/d approx target_value  => n approx target_value * d
        n = round(target_value * d)
        
        # Ensure numerator is within the bit limit and non-zero
        if 1 <= n <= limit:
            error = abs(n / d - target_value)
            if error < min_error:
                min_error = error
                best_fraction = (n, d)
                
    return best_fraction

def solve_titan_mass_problem():
    """
    Solves the Titan mass calculation problem step-by-step.
    """
    # 1. Define constants and given values as 5-bit friendly fractions
    # Max integer value for 5-bit unsigned is 2^5 - 1 = 31
    LIMIT = 31
    
    frac_4_3 = (4, 3)
    # The best 5-bit approximation for pi is 22/7
    frac_pi = (22, 7)
    # r = 0.5 cm -> 1/2. r^3 = (1/2)^3 = 1/8
    frac_r_cubed = (1, 8)
    # density = 0.9 kg/cm^3 -> 9/10
    frac_rho = (9, 10)

    print("Step 1: Set up the mass calculation with fractional components.")
    print(f"Mass = (4/3) * pi * r^3 * rho")
    print(f"       = ({frac_4_3[0]}/{frac_4_3[1]}) * ({frac_pi[0]}/{frac_pi[1]}) * ({frac_r_cubed[0]}/{frac_r_cubed[1]}) * ({frac_rho[0]}/{frac_rho[1]})")
    print("-" * 30)

    # 2. Perform the ideal calculation by multiplying all numerators and denominators
    # This theoretical step assumes we can handle larger intermediate numbers before final simplification.
    num = frac_4_3[0] * frac_pi[0] * frac_r_cubed[0] * frac_rho[0]
    den = frac_4_3[1] * frac_pi[1] * frac_r_cubed[1] * frac_rho[1]
    
    # We can simplify this ideal fraction using the greatest common divisor (GCD)
    common_divisor = math.gcd(num, den)
    ideal_num = num // common_divisor
    ideal_den = den // common_divisor

    print(f"Step 2: Perform the ideal calculation.")
    print(f"Mass = ({num}) / ({den})")
    print(f"Simplified, the ideal result is {ideal_num}/{ideal_den}.")
    print("-" * 30)
    
    # 3. Check if the ideal result is representable
    print(f"Step 3: Check if the result is valid on Titan.")
    if ideal_num > LIMIT or ideal_den > LIMIT:
        print(f"The result {ideal_num}/{ideal_den} is not representable, as {ideal_num} or {ideal_den} exceeds the 5-bit limit of {LIMIT}.")
    else:
        print(f"The result {ideal_num}/{ideal_den} is representable.")

    print("-" * 30)

    # 4. Find the best representable approximation
    target_decimal_value = ideal_num / ideal_den
    calculated_frac = find_best_approximation(target_decimal_value, LIMIT)
    
    print("Step 4: Find the closest representable fraction.")
    print(f"The decimal value of {ideal_num}/{ideal_den} is approx {target_decimal_value:.6f}.")
    print(f"The closest representable fraction (n,d <= {LIMIT}) is {calculated_frac[0]}/{calculated_frac[1]}.")
    calculated_mass = calculated_frac[0] / calculated_frac[1]
    print(f"This will be our calculated mass: {calculated_mass:.6f} kg.")
    print("-" * 30)

    # 5. Calculate the absolute error
    # True mass = (4/3) * math.pi * (0.5**3) * 0.9 = 0.15 * math.pi
    true_mass = 0.15 * math.pi
    absolute_error = abs(calculated_mass - true_mass)
    
    print("Step 5: Calculate the absolute error.")
    print(f"True mass (using high-precision pi) is approx {true_mass:.6f} kg.")
    print(f"The absolute error is |{calculated_mass:.6f} - {true_mass:.6f}| = {absolute_error:.6f}.")
    print("-" * 30)
    
    # 6. Final answer: error rounded to 0.001
    final_error = round(absolute_error, 3)
    print(f"The final error, rounded to 0.001, is {final_error}.")
    return final_error

# Run the calculation and print the final answer in the required format
final_answer = solve_titan_mass_problem()
print(f"<<<{final_answer}>>>")
