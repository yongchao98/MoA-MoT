import math

def is_prime(n):
    """Checks if a number is prime."""
    if n <= 1:
        return False
    for i in range(2, int(math.sqrt(n)) + 1):
        if n % i == 0:
            return False
    return True

def check_z_condition(r1, r2, r3):
    """
    Checks if z({r1, r2, r3}) == r2.
    The function z(C) returns the element x in C that maximizes |μ_C - μ_{C\\{x}}|.
    This is equivalent to maximizing |x - μ_C|.
    """
    c = [r1, r2, r3]
    mean = sum(c) / len(c)
    
    dist_r1 = abs(r1 - mean)
    dist_r2 = abs(r2 - mean)
    dist_r3 = abs(r3 - mean)
    
    # Return True if r2 has the maximum or equal-to-maximum distance from the mean.
    return dist_r2 >= dist_r1 and dist_r2 >= dist_r3

def solve():
    """
    Solves the problem by finding valid resistor values and calculating the max current.
    """
    R2 = 6
    found_resistors = []

    # From analysis, 6 must be the minimum resistor, and R3 > 10.
    # We search for the smallest prime R3 > 10.
    r3_candidate = 11
    while not found_resistors:
        if is_prime(r3_candidate):
            R3 = r3_candidate
            
            # From analysis, R1 is constrained by (R3+6)/2 < R1 < R3-2.
            min_r1 = math.floor((R3 + 6) / 2) + 1
            max_r1 = R3 - 3 # range is exclusive on the upper bound, R1 < R3-2 means up to R3-3

            for R1 in range(min_r1, max_r1 + 1):
                # Ensure resistors are distinct
                if R1 == R2 or R1 == R3:
                    continue

                # Check the z-condition
                if check_z_condition(R1, R2, R3):
                    # Found a valid set of resistors. We assume the first set found (with the
                    # smallest R values) will yield the maximum current.
                    found_resistors = [R1, R2, R3]
                    break
        if found_resistors:
            break
        r3_candidate += 1

    if not found_resistors:
        print("No solution found that satisfies all conditions.")
        return

    R1, R2, R3 = found_resistors
    print("Step 1 & 2: Finding Resistor Values")
    print("-" * 35)
    print(f"Given conditions lead to the smallest valid set of distinct integer resistors:")
    print(f"R1 = {R1} Ω")
    print(f"R2 = {R2} Ω")
    print(f"R3 = {R3} Ω (prime)")
    
    mean = (R1 + R2 + R3) / 3
    print("\nVerifying conditions:")
    print(f"1. z({R1}, {R2}, {R3}) = 6: The element furthest from the mean ({mean:.2f}) must be 6.")
    print(f"   |{R1} - {mean:.2f}| = {abs(R1-mean):.2f}")
    print(f"   |{R2} - {mean:.2f}| = {abs(R2-mean):.2f}")
    print(f"   |{R3} - {mean:.2f}| = {abs(R3-mean):.2f}")
    print(f"   Condition is met since {abs(R2-mean):.2f} is the maximum distance.")
    
    print(f"2. R3 - R1 > 2: {R3} - {R1} = {R3 - R1}, which is greater than 2. Condition is met.")

    print("\nStep 3 & 4: Analyzing Circuits and Calculating Max Current")
    print("-" * 60)
    print("Two circuit configurations can cause voltage across R3 to rise when R2 fails.")
    
    # Configuration 1: All resistors in parallel
    print("\nConfiguration A: All resistors in parallel")
    print("If R2 fails (open), the source current I_s is split between R1 and R3.")
    print(f"The voltage across R3 becomes 26 V, so I_s * (R1 || R3) = 26 V.")
    # I_s = 26 * (1/R1 + 1/R3)
    I_s_1_num = 26 * (R1 + R3)
    I_s_1_den = R1 * R3
    I_s_1 = I_s_1_num / I_s_1_den
    
    # When R2 is intact, this source current I_s is split between R1, R2, and R3.
    # The voltage across them is V_intact = I_s / (1/R1 + 1/R2 + 1/R3)
    R_eq_p_inv = 1/R1 + 1/R2 + 1/R3
    V_intact_1 = I_s_1 / R_eq_p_inv
    I3_1 = V_intact_1 / R3
    
    print("When R2 is intact, the current through R3 is:")
    print(f"I3_A = (I_s) / (R3 * (1/R1 + 1/R2 + 1/R3))")
    # Equation with numbers
    num = 156 * (R1 + R3)
    den = R3 * (R2*R3 + R1*R3 + R1*R2)
    print(f"I3_A = (156 * ({R1} + {R3})) / ({R3} * ({R2}*{R3} + {R1}*{R3} + {R1}*{R2}))")
    print(f"I3_A = {num} / {den} = 69 / 67 A ≈ {I3_1:.4f} A")

    # Configuration 2: R1 in series with R2 || R3
    print("\nConfiguration B: R1 in series with (R2 || R3)")
    print("If R2 fails (open), the source current I_s flows through R1 and R3 in series.")
    print(f"The voltage across R3 becomes 26 V, so I_s * R3 = 26 V.")
    I_s_2 = 26 / R3

    # When R2 is intact, I_s flows through R1 and then splits between R2 and R3.
    # We find the current through R3 using the current divider rule on the parallel part.
    # I3 = I_s * R2 / (R2 + R3)
    I3_2 = I_s_2 * R2 / (R2 + R3)
    
    print("When R2 is intact, the current through R3 is:")
    print(f"I3_B = I_s * R2 / (R2 + R3) = (26/R3) * R2 / (R2 + R3)")
    # Equation with numbers
    num = 26 * R2
    den = R3 * (R2 + R3)
    print(f"I3_B = (26 * {R2}) / ({R3} * ({R2} + {R3}))")
    print(f"I3_B = {num} / {den} = 12 / 19 A ≈ {I3_2:.4f} A")

    max_current = max(I3_1, I3_2)
    print("\n--- Conclusion ---")
    print(f"The two possible currents for R3 are {69/67:.4f} A and {12/19:.4f} A.")
    print(f"The maximum possible current through R3 is the larger of the two values.")
    print(f"Maximum Current = {max_current:.4f} A")

    return max_current

# Execute the solver function
max_current_value = solve()
# Format the final answer as requested
print(f"\nFinal Answer: The maximum possible current is 69/67 amperes.")
final_answer_float = 69 / 67
print(f"<<<{final_answer_float:.5f}>>>")
