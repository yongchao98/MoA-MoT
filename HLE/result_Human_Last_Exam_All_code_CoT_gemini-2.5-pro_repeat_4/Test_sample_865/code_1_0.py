import math

def format_equation_side(angles, t_index, side_indices):
    """Helper to format one side of the Kawasaki equation string."""
    parts = []
    for i in side_indices:
        if i == t_index:
            parts.append('t')
        else:
            # Use int for cleaner display if it's a whole number
            val = angles[i]
            parts.append(str(int(val)) if val == int(val) else str(val))
    return " + ".join(parts)

def solve_crease_pattern(pattern_data):
    """
    Analyzes a single-vertex crease pattern to find the angle 't' for flat-foldability.
    Prints the step-by-step analysis and returns the result.
    """
    pattern_num, pattern = pattern_data
    print(f"--- Analyzing Case {pattern_num} ---")
    
    # 1. Parse the input pattern into angles and folds
    angles = []
    folds = []
    t_index = -1
    for i, item in enumerate(pattern):
        if i % 2 == 0:  # Angle
            if item == 't':
                angles.append('t')
                t_index = len(angles) - 1
            else:
                angles.append(float(item))
        else:  # Fold type
            folds.append(item)

    # 2. Check Maekawa's Theorem
    num_m = folds.count('M')
    num_v = folds.count('V')
    print(f"Found {num_m} Mountain folds and {num_v} Valley folds.")
    if abs(num_m - num_v) != 2:
        print(f"Maekawa's Theorem Fails: The difference |{num_m} - {num_v}| = {abs(num_m - num_v)} is not 2.")
        return "none"
    print(f"Maekawa's Theorem Holds: |{num_m} - {num_v}| = {abs(num_m - num_v)} = 2.")

    # Handle cases with no variable angle 't'
    if t_index == -1:
        print("The pattern has no unknown angle 't' to solve for.")
        return "none"
        
    # 3. Solve for 't' using the sum-to-360 rule
    known_angles_sum = sum(a for a in angles if a != 't')
    t_from_sum = 360.0 - known_angles_sum
    
    angle_terms_str = " + ".join([str(int(a)) if a != 't' and a == int(a) else str(a) for a in angles])
    sum_eq_str = f"{angle_terms_str} = 360"
    print("\nApplying the rule that all angles must sum to 360 degrees:")
    print(f"  Equation: {sum_eq_str}")
    print(f"  Solving for t: {known_angles_sum} + t = 360  =>  t = {t_from_sum}")

    # 4. Solve for 't' using Kawasaki's Theorem
    odd_indices = list(range(0, len(angles), 2))
    even_indices = list(range(1, len(angles), 2))
    
    sum_odd_known = sum(angles[i] for i in odd_indices if i != t_index)
    sum_even_known = sum(angles[i] for i in even_indices if i != t_index)
    
    odd_side_str = format_equation_side(angles, t_index, odd_indices)
    even_side_str = format_equation_side(angles, t_index, even_indices)
    kawasaki_eq_str = f"{odd_side_str} = {even_side_str}"
    
    print("\nApplying Kawasaki's Theorem (alternating angle sums must be equal):")
    print(f"  Equation: {kawasaki_eq_str}")
    
    if t_index in odd_indices:
        t_from_kawasaki = sum_even_known - sum_odd_known
        print(f"  Solving for t: {sum_odd_known} + t = {sum_even_known}  =>  t = {t_from_kawasaki}")
    else: # t_index is in even_indices
        t_from_kawasaki = sum_odd_known - sum_even_known
        print(f"  Solving for t: {sum_odd_known} = {sum_even_known} + t  =>  t = {t_from_kawasaki}")

    # 5. Compare results and determine the final answer
    print("\nConclusion:")
    if t_from_sum > 0 and math.isclose(t_from_sum, t_from_kawasaki):
        result = int(round(t_from_sum))
        print(f"The two methods agree on a valid angle t = {result}.")
        return result
    else:
        print(f"The required values for t are inconsistent ({t_from_sum} vs {t_from_kawasaki}).")
        print("Therefore, no solution exists that makes the pattern flat-foldable.")
        return "none"

if __name__ == '__main__':
    patterns_to_solve = [
        (1, [100, 'M', 62, 'V', 22, 'M', 't', 'V', 33, 'M', 90, 'V']),
        (2, [90, 'M', 120, 'M', 60, 'M', 90, 'M']),
        (3, [60, 'V', 60, 'M', 120, 'M', 't', 'M']),
        (4, [77, 'M', 15, 'M', 50, 'V', 33, 'M', 't', 'V', 130, 'M'])
    ]
    
    final_results = []
    for p_data in patterns_to_solve:
        result = solve_crease_pattern(p_data)
        final_results.append(result)
        print("-" * 30)
        
    print("\nFinal comma-separated list of values for t:")
    print(f"[{','.join(map(str, final_results))}]")