def solve_crease_pattern(pattern_data, pattern_id):
    """
    Analyzes a single-vertex crease pattern for flat-foldability,
    solves for the unknown angle 't', and prints the reasoning.
    """
    print(f"--- Analyzing Pattern {pattern_id} ---")
    print(f"Pattern: {pattern_data}")

    angles = [p for i, p in enumerate(pattern_data) if i % 2 == 0]
    folds = [p for i, p in enumerate(pattern_data) if i % 2 == 1]

    # 1. Maekawa's Theorem Check: |#M - #V| = 2
    num_m = folds.count('M')
    num_v = folds.count('V')
    if abs(num_m - num_v) != 2:
        print(f"Status: Fails Maekawa's Theorem.")
        print(f"Reason: The number of mountain folds (#M={num_m}) and valley folds (#V={num_v}) must differ by exactly 2.")
        print(f"Equation: |{num_m} - {num_v}| = {abs(num_m - num_v)} != 2")
        print(f"Result for Pattern {pattern_id}: none")
        return 'none'

    print("Status: Maekawa's Theorem holds.")

    # 2. Kawasaki's Theorem Check: Sum of alternating angles is 180.
    odd_angles_list = []
    even_angles_list = []
    t_is_in_odd_group = False
    
    # Handle patterns without 't' (like pattern 2) by treating it as a check
    has_t = 't' in angles

    for i, angle in enumerate(angles):
        is_odd_angle = (i % 2 == 0) # 1st, 3rd, 5th... angles are at even indices 0, 2, 4...
        if is_odd_angle:
            odd_angles_list.append(angle)
            if angle == 't':
                t_is_in_odd_group = True
        else:
            even_angles_list.append(angle)

    group_to_check = even_angles_list if t_is_in_odd_group else odd_angles_list
    group_to_solve = odd_angles_list if t_is_in_odd_group else even_angles_list
    
    # Check the sum of the group of angles that does not contain 't'
    known_angles_sum = sum([a for a in group_to_check if isinstance(a, int) or isinstance(a, float)])
    
    if known_angles_sum != 180:
        equation_str = " + ".join(map(str, group_to_check))
        print("Status: Fails Kawasaki's Theorem.")
        print("Reason: The sum of alternating angles without 't' must equal 180.")
        print(f"Equation: {equation_str} = {known_angles_sum} != 180")
        print(f"Result for Pattern {pattern_id}: none")
        return 'none'

    # If the check passes, solve for 't' using the other group of angles
    if not has_t: # Should not happen based on prompt, but handles patterns without 't'
        print(f"Status: Fails Kawasaki's Theorem (no 't' to solve for).")
        print(f"Result for Pattern {pattern_id}: none")
        return 'none'

    known_in_solve_group = [a for a in group_to_solve if a != 't']
    sum_known_solve = sum(known_in_solve_group)
    t_value = 180 - sum_known_solve

    equation_parts = [str(a) for a in known_in_solve_group] + ['t']
    equation_str = " + ".join(equation_parts)

    print("Status: Kawasaki's Theorem holds.")
    print(f"Reason: Solving for 't' using the alternating angles.")
    print(f"Equation: {equation_str} = 180")
    
    if t_value <= 0:
        print(f"Result for Pattern {pattern_id}: none (solving gives non-positive angle {t_value})")
        return 'none'

    print(f"Result for Pattern {pattern_id}: {t_value}")
    return t_value

if __name__ == '__main__':
    patterns = [
        [100, 'M', 62, 'V', 22, 'M', 't', 'V', 33, 'M', 90, 'V'],
        [90, 'M', 120, 'M', 60, 'M', 90, 'M'],
        [60, 'V', 60, 'M', 120, 'M', 't', 'M'],
        [77, 'M', 15, 'M', 50, 'V', 33, 'M', 't', 'V', 130, 'M']
    ]

    final_results = []
    for i, p in enumerate(patterns):
        result = solve_crease_pattern(p, i + 1)
        final_results.append(result)
        print("-" * 25)

    print("\nFinal comma-separated list of values:")
    # Format the list with 'none' as a string for clarity
    formatted_results = [res if res != 'none' else 'none' for res in final_results]
    print(formatted_results)