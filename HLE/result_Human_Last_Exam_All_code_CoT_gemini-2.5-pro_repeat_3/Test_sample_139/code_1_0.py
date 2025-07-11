import math

def is_prime(n):
    """Checks if a number is prime."""
    if n <= 1:
        return False
    for i in range(2, int(math.sqrt(n)) + 1):
        if n % i == 0:
            return False
    return True

def find_max_current():
    """
    Finds the resistor values R1 and R3 that satisfy the problem's conditions
    and maximize the current through R3 in the intact circuit.
    """
    max_i3 = 0.0
    best_r1 = None
    best_r3 = None
    r2 = 6

    # The current function decreases as R1 and R3 increase, so we only need to search a reasonable range.
    for r1_candidate in range(1, 100):
        for r3_candidate in range(1, 100):
            # Condition 1: R1, R2, R3 are distinct integers
            if r1_candidate == r2 or r3_candidate == r2 or r1_candidate == r3_candidate:
                continue

            # Condition 2: R3 is a prime number
            if not is_prime(r3_candidate):
                continue

            # Condition 3: R3 - R1 > 2
            if not (r3_candidate - r1_candidate > 2):
                continue
                
            # Condition 4: z(R1, 6, R3) = 6, which implies R1 > 6, R3 > 6 and two inequalities
            # These two inequalities were derived from the definition of z(C)
            # R3 + 6 <= 2*R1  and  R1 + 6 <= 2*R3
            if not (r3_candidate + r2 <= 2 * r1_candidate and r1_candidate + r2 <= 2 * r3_candidate):
                continue

            # If all conditions are met, calculate the current
            # I_R3 = [26 * (R1 + R3) * R2] / [R3 * (R1*R2 + R1*R3 + R2*R3)]
            numerator = 26 * (r1_candidate + r3_candidate) * r2
            denominator = r3_candidate * (r1_candidate * r2 + r1_candidate * r3_candidate + r2 * r3_candidate)
            
            if denominator == 0:
                continue

            current_i3 = numerator / denominator

            if current_i3 > max_i3:
                max_i3 = current_i3
                best_r1 = r1_candidate
                best_r3 = r3_candidate

    return best_r1, best_r3, r2, max_i3

# Find the optimal values
R1, R3, R2, max_current = find_max_current()

if R1 is not None:
    print(f"The search found the optimal resistor values that satisfy all conditions:")
    print(f"R1 = {R1} ohms")
    print(f"R2 = {R2} ohms")
    print(f"R3 = {R3} ohms (prime)")
    print("\nThese values result in the maximum possible current through R3.")
    
    print("\nThe equation for the current through R3 when the circuit is intact is:")
    print("I_R3 = [26 * (R1 + R3) * R2] / [R3 * (R1*R2 + R1*R3 + R2*R3)]")
    
    print("\nPlugging in the optimal values:")
    # Show the calculation step-by-step as requested
    num_calc_1 = f"26 * ({R1} + {R3}) * {R2}"
    num_val_1 = 26 * (R1 + R3) * R2
    
    den_calc_1 = f"{R3} * ({R1}*{R2} + {R1}*{R3} + {R2}*{R3})"
    
    den_val_part1 = R1 * R2
    den_val_part2 = R1 * R3
    den_val_part3 = R2 * R3
    den_calc_2 = f"{R3} * ({den_val_part1} + {den_val_part2} + {den_val_part3})"
    den_val_2 = R3 * (den_val_part1 + den_val_part2 + den_val_part3)
    
    
    print(f"I_R3 = [{num_calc_1}] / [{den_calc_1}]")
    print(f"I_R3 = [{num_val_1}] / [{den_calc_2}]")
    print(f"I_R3 = {num_val_1} / {den_val_2}")

    # Simplify fraction using GCD
    common_divisor = math.gcd(num_val_1, den_val_2)
    simplified_num = num_val_1 // common_divisor
    simplified_den = den_val_2 // common_divisor
    
    print(f"I_R3 = {simplified_num} / {simplified_den} A")
    
    print(f"\nThe maximum possible current through R3 is approximately {max_current:.4f} A.")

else:
    print("No solution found within the search range.")

print(f"\n<<<{max_current:.4f}>>>")