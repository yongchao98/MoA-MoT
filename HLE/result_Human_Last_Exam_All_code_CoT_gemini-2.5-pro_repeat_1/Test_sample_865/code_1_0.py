def solve_flat_foldability():
    """
    Analyzes single-vertex crease patterns to find the angle 't' that allows for flat-folding.
    It checks Maekawa's and Kawasaki's theorems.
    """
    patterns = [
        [100, 'M', 62, 'V', 22, 'M', 't', 'V', 33, 'M', 90, 'V'],
        [90, 'M', 120, 'M', 60, 'M', 90, 'M'],
        [60, 'V', 60, 'M', 120, 'M', 't', 'M'],
        [77, 'M', 15, 'M', 50, 'V', 33, 'M', 't', 'V', 130, 'M']
    ]

    final_results = []

    for i, pattern_data in enumerate(patterns):
        angles = []
        folds = []
        t_index = -1
        known_angle_sum = 0.0

        # Parse the pattern into separate lists for angles and folds
        for j in range(0, len(pattern_data), 2):
            angle_val = pattern_data[j]
            if isinstance(angle_val, str) and angle_val == 't':
                t_index = len(angles)
                angles.append('t')
            else:
                angles.append(float(angle_val))
                known_angle_sum += float(angle_val)
        
        for j in range(1, len(pattern_data), 2):
            folds.append(pattern_data[j])

        # 1. Maekawa's Theorem: |#M - #V| = 2
        num_m = folds.count('M')
        num_v = folds.count('V')
        if abs(num_m - num_v) != 2:
            final_results.append('none')
            continue

        # If no 't' is present, we cannot solve for it as per the prompt.
        if t_index == -1:
            final_results.append('none')
            continue

        # 2. Kawasaki's Theorem
        # From Sum of Angles: sum(all) = 360
        t_sum = 360.0 - known_angle_sum
        
        # From Alternating Sum: a1 - a2 + a3 - ... = 0
        alt_sum_known = 0.0
        for j, angle in enumerate(angles):
            if j != t_index:
                alt_sum_known += angle * ((-1)**j)
        
        # Deriving t from the alternating sum equation: t * (-1)^t_index + alt_sum_known = 0
        t_kawasaki = -alt_sum_known * ((-1)**t_index)

        # 3. Check for consistency and validity
        # The angle must be positive and both theorems must yield the same value for t.
        # A small tolerance is used for floating-point comparison.
        if t_sum <= 1e-9 or t_kawasaki <= 1e-9 or abs(t_sum - t_kawasaki) > 1e-9:
            final_results.append('none')
        else:
            solution_t = int(round(t_sum))
            final_results.append(solution_t)
            
            # As requested, output the equations for the valid solution
            final_angles = list(angles)
            final_angles[t_index] = solution_t

            print(f"Solution found for pattern {i+1} with t = {solution_t}")

            # Print the sum of angles equation
            sum_eq_parts = [str(int(a)) for a in final_angles]
            print(f"Sum of angles verification: {' + '.join(sum_eq_parts)} = {int(sum(final_angles))}")
            
            # Print the alternating sum of angles equation
            alt_sum_eq_parts = []
            calculated_alt_sum = 0
            for j, angle in enumerate(final_angles):
                term = int(angle)
                if j > 0:
                    if j % 2 != 0:
                        alt_sum_eq_parts.append(f"- {term}")
                    else:
                        alt_sum_eq_parts.append(f"+ {term}")
                else:
                    alt_sum_eq_parts.append(f"{term}")
                calculated_alt_sum += term * ((-1)**j)
            
            print(f"Alternating sum verification: {' '.join(alt_sum_eq_parts)} = {calculated_alt_sum}")
            print("-" * 20)

    # Format the final list for the final answer
    result_str = ','.join([str(r) for r in final_results])
    print(f"\nFinal Result List: [{result_str}]")


solve_flat_foldability()