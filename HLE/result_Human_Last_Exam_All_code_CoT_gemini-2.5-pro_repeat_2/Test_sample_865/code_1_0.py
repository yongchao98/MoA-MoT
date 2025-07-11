def solve_crease_pattern(pattern_data, case_num):
    """
    Analyzes a single-vertex crease pattern to find the value of angle 't'
    that makes it flat-foldable.
    """
    print(f"\n--- Analyzing Case {case_num} ---")
    print(f"Pattern: {pattern_data}")

    # 1. Parse pattern into angles and folds
    angles = []
    folds = []
    t_index = -1
    has_t = 't' in pattern_data

    for i, item in enumerate(pattern_data):
        if i % 2 == 0:  # Angle
            if item == 't':
                angles.append(None)
                t_index = len(angles) - 1
            else:
                angles.append(item)
        else:  # Fold
            folds.append(item)

    # If no 't' is present, we cannot solve for a value.
    if not has_t:
        print("Pattern does not contain an unknown angle 't' to solve for.")
        print("Result: none")
        return 'none'
        
    # 2. Check Maekawa's Theorem
    num_m = folds.count('M')
    num_v = folds.count('V')
    print(f"Maekawa's theorem: |#M - #V| = |{num_m} - {num_v}| = {abs(num_m - num_v)}")
    if abs(num_m - num_v) != 2:
        print("Condition not met (difference is not 2).")
        print("Result: none")
        return 'none'
    print("Condition met.")

    # 3. Check Kawasaki's Theorem by solving for 't'
    print("Kawasaki's theorem:")

    # 3a. Solve for 't' using Sum of Angles = 360
    known_angles = [a for a in angles if a is not None]
    sum_known_angles = sum(known_angles)
    t_sum = 360 - sum_known_angles

    equation_parts_sum = [str(a) for a in known_angles] + ['t']
    print(f"  1. Sum of angles: {' + '.join(equation_parts_sum)} = 360")
    print(f"     => {sum_known_angles} + t = 360")
    print(f"     => t = {t_sum}")

    # 3b. Solve for 't' using Alternating Sum = 0
    sum_odd_placed = 0
    sum_even_placed = 0
    odd_parts = []
    even_parts = []

    # 0-indexed: index 0 is 1st (odd), index 1 is 2nd (even), etc.
    t_is_odd_placed = (t_index % 2 == 0)

    for i, angle in enumerate(angles):
        is_odd_placed = (i % 2 == 0)
        if angle is None: # This is 't'
            if is_odd_placed:
                odd_parts.append('t')
            else:
                even_parts.append('t')
            continue

        if is_odd_placed:
            sum_odd_placed += angle
            odd_parts.append(str(angle))
        else:
            sum_even_placed += angle
            even_parts.append(str(angle))

    print(f"  2. Alternating sum: {' + '.join(odd_parts)} = {' + '.join(even_parts)}")
    if t_is_odd_placed:
        t_alt = sum_even_placed - sum_odd_placed
        print(f"     => {sum_odd_placed} + t = {sum_even_placed}")
        print(f"     => t = {t_alt}")
    else:
        t_alt = sum_odd_placed - sum_even_placed
        print(f"     => {sum_odd_placed} = {sum_even_placed} + t")
        print(f"     => t = {t_alt}")

    # 4. Compare results and conclude
    print("\nConclusion:")
    if t_sum != t_alt:
        print(f"Values for t are inconsistent ({t_sum} vs {t_alt}). The pattern cannot be flat-folded.")
        print("Result: none")
        return 'none'
    
    if t_sum <= 0:
        print(f"Resulting angle t={t_sum} must be positive.")
        print("Result: none")
        return 'none'
        
    print(f"Both equations give t = {t_sum}. All conditions are met.")
    print(f"Result: {t_sum}")
    return t_sum

def main():
    """
    Main function to solve the four crease pattern problems.
    """
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

    # Format the final list as a string like "[72,33,none,5]"
    result_strings = [str(r) if r != 'none' else 'none' for r in results]
    formatted_results = f"[{','.join(result_strings)}]"
    
    print("\n" + "="*30)
    print("Final Answer Summary")
    print("="*30)
    print(formatted_results)
    print(f"<<<{formatted_results}>>>")

if __name__ == '__main__':
    main()