import math

def solve_crease_pattern(pattern_data, pattern_id):
    """
    Analyzes a single-vertex crease pattern to find the angle 't' that
    makes it flat-foldable. It prints its reasoning.

    Args:
        pattern_data: A list representing the crease pattern.
        pattern_id: An integer for identifying the pattern in the output.

    Returns:
        The numeric value of 't' if a solution exists, otherwise the string 'none'.
    """
    print(f"--- Analyzing Pattern {pattern_id} ---")
    
    angles_in = []
    folds = []
    t_index = -1
    known_angle_sum = 0
    has_t = False

    # 1. Parse the input pattern
    for i in range(0, len(pattern_data), 2):
        angle = pattern_data[i]
        fold = pattern_data[i+1]
        folds.append(fold)
        if isinstance(angle, str) and angle == 't':
            has_t = True
            t_index = i // 2
            angles_in.append('t')
        else:
            angles_in.append(angle)
            known_angle_sum += angle
            
    # Check for even number of creases
    if len(angles_in) % 2 != 0:
        print("Failed: The number of creases must be even.")
        return 'none'

    # 2. Check Maekawa's Theorem: |N_M - N_V| = 2
    n_m = folds.count('M')
    n_v = folds.count('V')
    if abs(n_m - n_v) != 2:
        print(f"Failed Maekawa's Theorem: The difference between mountain and valley folds is not 2.")
        print(f"Equation: |{n_m} (M) - {n_v} (V)| = {abs(n_m - n_v)}")
        return 'none'
    print(f"Passed Maekawa's Theorem: |{n_m} (M) - {n_v} (V)| = 2")

    # 3. Solve for t by summing angles to 360
    t_value = None
    if has_t:
        t_value = 360 - known_angle_sum
        known_angles_str = " + ".join(map(str, [a for a in angles_in if a != 't']))
        print(f"Solving for t: {known_angles_str} + t = 360")
        print(f"Result: t = {t_value}")
        if t_value <= 0:
            print(f"Failed: Angle 't' must be positive, but calculated t = {t_value}.")
            return 'none'
        angles = [t_value if a == 't' else a for a in angles_in]
    else:
        if not math.isclose(known_angle_sum, 360):
            print(f"Failed: Sum of angles is {known_angle_sum}, not 360, and no 't' to adjust.")
            return 'none'
        angles = angles_in
    
    # 4. Check Kawasaki's Theorem: Sum of alternating angles are equal
    sum_odd_indices = sum(angles[i] for i in range(0, len(angles), 2))
    sum_even_indices = sum(angles[i] for i in range(1, len(angles), 2))
    
    odd_angles_str = " + ".join(map(str, [angles[i] for i in range(0, len(angles), 2)]))
    even_angles_str = " + ".join(map(str, [angles[i] for i in range(1, len(angles), 2)]))

    if not math.isclose(sum_odd_indices, sum_even_indices):
        print(f"Failed Kawasaki's Theorem: Alternating angle sums are not equal.")
        print(f"Equation: {odd_angles_str} = {sum_odd_indices} != {even_angles_str} = {sum_even_indices}")
        return 'none'
    
    print("Passed Kawasaki's Theorem: Alternating angle sums are equal.")
    print(f"Equation: {odd_angles_str} = {sum_odd_indices} and {even_angles_str} = {sum_even_indices}")

    # 5. Conclusion
    if has_t:
        print(f"Success! A valid angle t = {t_value} exists.")
        return int(t_value) if t_value.is_integer() else t_value
    else:
        # This case applies if a pattern without 't' is checked.
        print("The pattern is flat-foldable, but no angle 't' was specified to solve for.")
        return 'none'


if __name__ == '__main__':
    patterns = [
        [100, 'M', 62, 'V', 22, 'M', 't', 'V', 33, 'M', 90, 'V'],
        [90, 'M', 120, 'M', 60, 'M', 90, 'M'],
        [60, 'V', 60, 'M', 120, 'M', 't', 'M'],
        [77, 'M', 15, 'M', 50, 'V', 33, 'M', 't', 'V', 130, 'M']
    ]

    results = []
    for i, p in enumerate(patterns):
        result = solve_crease_pattern(p, i + 1)
        results.append(result)
        print("") # Add a newline for readability

    # Format the final output string
    final_answer_str = f"[{','.join(str(r) for r in results)}]"
    print(f"Final comma-separated list: {final_answer_str}")
    print(f"<<<{final_answer_str}>>>")