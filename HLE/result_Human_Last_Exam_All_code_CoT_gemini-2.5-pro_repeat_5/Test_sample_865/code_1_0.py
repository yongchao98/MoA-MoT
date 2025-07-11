import math

def solve_flat_fold_angle(patterns):
    """
    Analyzes four single-vertex crease patterns to find the angle 't' that
    makes them flat-foldable.

    It prints a step-by-step analysis for each pattern and a final summary list.
    """
    final_results = []

    for i, pattern_data in enumerate(patterns):
        case_num = i + 1
        print(f"\n--- Analyzing Case {case_num}: {pattern_data} ---")

        angles = []
        folds = []
        t_info = {'present': False, 'index': -1}

        # 1. Parse the input pattern into angles and folds
        for j, item in enumerate(pattern_data):
            if j % 2 == 0:  # It's an angle
                if item == 't':
                    t_info['present'] = True
                    t_info['index'] = len(angles)
                    angles.append('t')
                else:
                    angles.append(float(item))
            else:  # It's a fold type
                folds.append(item)

        # 2. Check Maekawa's Theorem: |M - V| = 2
        num_m = folds.count('M')
        num_v = folds.count('V')

        if abs(num_m - num_v) != 2:
            print(f"Condition 1: Maekawa's Theorem (|M - V| = 2)")
            print(f"FAIL: The number of Mountain folds ({num_m}) and Valley folds ({num_v}) does not differ by 2.")
            print("Conclusion: The pattern is not flat-foldable.")
            final_results.append('none')
            continue

        print(f"Condition 1: Maekawa's Theorem (|M - V| = 2)")
        print(f"PASS: |{num_m} M - {num_v} V| = {abs(num_m - num_v)}.")

        # 3. Check Kawasaki's Theorem: Sum of alternating angles = 180°
        print(f"\nCondition 2: Kawasaki's Theorem (Alternating angle sums must be 180°)")
        sum_odd_angles = 0
        sum_even_angles = 0
        t_is_in_odd_group = False
        
        known_odd_angles = []
        known_even_angles = []

        for j, angle in enumerate(angles):
            is_odd_indexed = (j + 1) % 2 != 0
            if angle == 't':
                if is_odd_indexed:
                    t_is_in_odd_group = True
            else:
                if is_odd_indexed:
                    sum_odd_angles += angle
                    known_odd_angles.append(str(int(angle)))
                else:
                    sum_even_angles += angle
                    known_even_angles.append(str(int(angle)))

        # If 't' is not present, we can't solve for it.
        if not t_info['present']:
             print("FAIL: The pattern has no variable angle 't' to solve for.")
             print("Conclusion: No value can be specified.")
             final_results.append('none')
             continue

        # Solve for 't' using the Kawasaki condition
        solution_found = False
        t_value = None

        if t_is_in_odd_group:
            # The sum of even angles must be 180 for a solution to exist
            even_sum_str = " + ".join(known_even_angles)
            print(f"Checking the sum of the known even-indexed angles: {even_sum_str} = {sum_even_angles}")
            if not math.isclose(sum_even_angles, 180):
                print("FAIL: This sum is not 180, so the pattern cannot be flat-foldable.")
                final_results.append('none')
            else:
                print("PASS: The sum is 180.")
                t_value = 180 - sum_odd_angles
                odd_sum_str = " + ".join(known_odd_angles)
                print(f"Solving for 't' using the odd-indexed angles:")
                print(f"Equation: {odd_sum_str} + t = 180")
                print(f"{sum_odd_angles} + t = 180")
                print(f"t = 180 - {sum_odd_angles}")
                solution_found = True
        else: # t is in the even group
            # The sum of odd angles must be 180 for a solution to exist
            odd_sum_str = " + ".join(known_odd_angles)
            print(f"Checking the sum of the known odd-indexed angles: {odd_sum_str} = {sum_odd_angles}")
            if not math.isclose(sum_odd_angles, 180):
                print("FAIL: This sum is not 180, so the pattern cannot be flat-foldable.")
                final_results.append('none')
            else:
                print("PASS: The sum is 180.")
                t_value = 180 - sum_even_angles
                even_sum_str = " + ".join(known_even_angles)
                print(f"Solving for 't' using the even-indexed angles:")
                print(f"Equation: {even_sum_str} + t = 180")
                print(f"{sum_even_angles} + t = 180")
                print(f"t = 180 - {sum_even_angles}")
                solution_found = True

        if solution_found:
            # Angles must be positive
            if t_value > 0:
                print(f"Conclusion: The required angle is t = {t_value}")
                final_results.append(int(t_value) if t_value.is_integer() else t_value)
            else:
                print(f"FAIL: Solving for t results in a non-positive angle ({t_value}), which is not possible.")
                final_results.append('none')

    print("\n" + "="*40)
    print("Final Summary")
    print("="*40)
    print("The values for 't' for each case are:")
    print(final_results)


if __name__ == '__main__':
    # The four single-vertex crease patterns from the problem
    patterns_to_solve = [
        [100, 'M', 62, 'V', 22, 'M', 't', 'V', 33, 'M', 90, 'V'],
        [90, 'M', 120, 'M', 60, 'M', 90, 'M'],
        [60, 'V', 60, 'M', 120, 'M', 't', 'M'],
        [77, 'M', 15, 'M', 50, 'V', 33, 'M', 't', 'V', 130, 'M']
    ]
    solve_flat_fold_angle(patterns_to_solve)