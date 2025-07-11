import math

def solve_crease_pattern_assignments():
    """
    Determines the number of valid flat-foldable assignments for a given crease pattern.
    """
    pattern = [60, 'M', 30, '?', 50, '?', 70, 'V', 150, '?']

    # 1. Parse the input
    angles = [p for p in pattern if isinstance(p, (int, float))]
    creases = [p for p in pattern if isinstance(p, str)]
    num_creases = len(creases)
    
    print("Analyzing the given crease pattern for a single vertex.")
    print("-" * 60)

    # 2. Check fundamental conditions for flat-foldability.
    # We will check two primary geometric conditions. A failure in either means 0 valid assignments.
    
    # Condition A: Number of creases must be even.
    print("Condition 1: The number of creases must be an even number.")
    is_crease_count_valid = (num_creases % 2 == 0)
    if not is_crease_count_valid:
        print(f"-> FAILED: The pattern has {num_creases} creases, which is an odd number.")
    else:
        print(f"-> PASSED: The pattern has {num_creases} creases, which is an even number.")

    # Condition B: Kawasaki's Theorem.
    print("\nCondition 2: The sum of alternating angles must be 180 degrees.")
    
    # We form the two alternating sums of angles.
    sum1_angles = angles[0::2]
    sum2_angles = angles[1::2]
    sum1 = sum(sum1_angles)
    sum2 = sum(sum2_angles)

    # Print the equations as requested
    sum1_eq_str = " + ".join(map(str, sum1_angles)) + f" = {sum1}"
    print(f"Sum of one set of alternating angles: {sum1_eq_str}")
    
    sum2_eq_str = " + ".join(map(str, sum2_angles)) + f" = {sum2}"
    print(f"Sum of the other set of alternating angles: {sum2_eq_str}")

    is_kawasaki_valid = math.isclose(sum1, 180) and math.isclose(sum2, 180)
    if not is_kawasaki_valid:
        print(f"-> FAILED: The sums are not 180 degrees.")
    else:
        print(f"-> PASSED: The sums are 180 degrees.")

    print("-" * 60)
    
    # 3. Conclude based on the checks.
    if not is_crease_count_valid or not is_kawasaki_valid:
        print("\nConclusion: The fundamental geometric properties for a flat-foldable vertex are not met.")
        print("Therefore, no matter how the unassigned creases are assigned, the pattern cannot be folded flat.")
        final_answer = 0
    else:
        # This part is for a valid geometry, which is not the case here.
        # It would calculate combinations based on Maekawa's Theorem |M-V|=2.
        print("\nConclusion: The geometric properties are met. Now checking Maekawa's Theorem for crease assignments.")
        num_M_known = creases.count('M')
        num_V_known = creases.count('V')
        num_unknown = creases.count('?')
        
        count = 0
        from math import comb
        # Iterate through all possible assignments for '?'
        # i represents the number of '?' assigned as 'M'
        for i in range(num_unknown + 1):
            num_M_total = num_M_known + i
            num_V_total = num_V_known + (num_unknown - i)
            
            if abs(num_M_total - num_V_total) == 2:
                count += comb(num_unknown, i)
        final_answer = count
    
    print(f"\nThe total number of different flat-foldable assignments is {final_answer}.")

solve_crease_pattern_assignments()