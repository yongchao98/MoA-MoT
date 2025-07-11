import math

def solve_crease_pattern(pattern_data, pattern_num):
    """
    Analyzes a single-vertex crease pattern to find the angle 't' for flat-foldability.
    """
    print(f"--- Analyzing Pattern {pattern_num} ---")
    print(f"Pattern: {pattern_data}")

    angles = []
    folds = []
    t_angle_index = -1  # 0-based index of 't' in the angles list

    for i, item in enumerate(pattern_data):
        if i % 2 == 0:  # It's an angle
            if item == 't':
                t_angle_index = len(angles)
            angles.append(item)
        else:  # It's a fold direction
            folds.append(item)

    # 1. Check Maekawa's Theorem: |M - V| = 2
    num_m = folds.count('M')
    num_v = folds.count('V')
    print(f"Maekawa's Theorem Check: Number of Mountains (M) = {num_m}, Number of Valleys (V) = {num_v}")
    
    if abs(num_m - num_v) != 2:
        print(f"Result: |M - V| = {abs(num_m - num_v)}, which is not 2. The pattern cannot be flat-folded.")
        return "none"
    print("Result: |M - V| = 2. Maekawa's Theorem is satisfied.")

    # 2. Check Angle Conditions if 't' is present
    if t_angle_index == -1:
        # This case applies to pattern 2. Since it already failed Maekawa's,
        # we can definitively say 'none'. This part handles hypothetical cases.
        print("Result: No variable 't' to solve for. Cannot determine a value for t.")
        return "none"

    # 2a. Solve for 't' using the Angle Sum Condition (sum = 360)
    known_angles = [a for a in angles if isinstance(a, int)]
    known_sum = sum(known_angles)
    t_from_sum = 360 - known_sum
    
    equation_str_sum = " + ".join(map(str, known_angles)) + " + t = 360"
    print(f"\nAngle Sum Equation: {equation_str_sum}")
    print(f"Solving for t: t = 360 - {known_sum} => t = {t_from_sum}")

    # 2b. Solve for 't' using Kawasaki's Theorem
    odd_sum_known = 0
    even_sum_known = 0
    
    if t_angle_index % 2 == 0: # 't' is an odd-indexed angle (1st, 3rd, 5th...)
        t_is_in_odd_group = True
    else: # 't' is an even-indexed angle (2nd, 4th, 6th...)
        t_is_in_odd_group = False

    # Collect known angle sums
    odd_angles_known_list = [a for i, a in enumerate(angles) if i % 2 == 0 and i != t_angle_index]
    even_angles_known_list = [a for i, a in enumerate(angles) if i % 2 != 0 and i != t_angle_index]
    odd_sum_known = sum(odd_angles_known_list)
    even_sum_known = sum(even_angles_known_list)

    # Build and solve Kawasaki's equation
    odd_str = " + ".join(map(str, odd_angles_known_list))
    even_str = " + ".join(map(str, even_angles_known_list))
    
    print("\nKawasaki's Theorem Equation (Sum of odd angles = Sum of even angles):")
    if t_is_in_odd_group:
        t_from_kawasaki = even_sum_known - odd_sum_known
        final_odd_str = (odd_str + " + t") if odd_str else "t"
        print(f"Equation: {final_odd_str} = {even_str}")
        print(f"Solving for t: t = {even_sum_known} - {odd_sum_known} => t = {t_from_kawasaki}")
    else:
        t_from_kawasaki = odd_sum_known - even_sum_known
        final_even_str = (even_str + " + t") if even_str else "t"
        print(f"Equation: {odd_str} = {final_even_str}")
        print(f"Solving for t: t = {odd_sum_known} - {even_sum_known} => t = {t_from_kawasaki}")

    # 3. Final Conclusion
    print("\nConclusion:")
    if not math.isclose(t_from_sum, t_from_kawasaki):
        print(f"The two required values for t ({t_from_sum} and {t_from_kawasaki}) are not equal.")
        print("Result: No single value of t satisfies all conditions.")
        return "none"
    
    if t_from_sum <= 0:
        print(f"The calculated value t = {t_from_sum} is not a positive angle.")
        print("Result: No valid angle t exists.")
        return "none"
        
    print(f"Both methods yield a consistent, positive value: t = {t_from_sum}.")
    print("Result: The pattern is flat-foldable with this angle.")
    return int(t_from_sum)


# Define the four crease patterns from the problem
patterns = [
    [100, 'M', 62, 'V', 22, 'M', 't', 'V', 33, 'M', 90, 'V'],
    [90, 'M', 120, 'M', 60, 'M', 90, 'M'],
    [60, 'V', 60, 'M', 120, 'M', 't', 'M'],
    [77, 'M', 15, 'M', 50, 'V', 33, 'M', 't', 'V', 130, 'M']
]

# Process each pattern and collect the results
final_results = []
for i, p in enumerate(patterns):
    result = solve_crease_pattern(p, i + 1)
    final_results.append(result)
    print("\n" + "="*40 + "\n")

# Print the final comma-separated list
final_output_str = "[" + ",".join(map(str, final_results)) + "]"
print(f"Final comma-separated list of values: {final_output_str}")
<<<[none,none,120,none]>>>