def solve_flat_foldability():
    """
    Solves for the unknown angle 't' in four single-vertex crease patterns
    and prints the justification for each case. Finally, it prints the
    results in a comma-separated list.
    """
    patterns = [
        [100, 'M', 62, 'V', 22, 'M', 't', 'V', 33, 'M', 90, 'V'],
        [90, 'M', 120, 'M', 60, 'M', 90, 'M'],
        [60, 'V', 60, 'M', 120, 'M', 't', 'M'],
        [77, 'M', 15, 'M', 50, 'V', 33, 'M', 't', 'V', 130, 'M']
    ]

    results = []

    for i, pattern in enumerate(patterns):
        case_num = i + 1
        angles_with_t = []
        folds = []

        # Separate angles and folds from the pattern
        for j, item in enumerate(pattern):
            if j % 2 == 0:
                angles_with_t.append(item)
            else:
                folds.append(item)
        
        # Check 1: Maekawa's Theorem
        n_M = folds.count('M')
        n_V = folds.count('V')

        if abs(n_M - n_V) != 2:
            print(f"For pattern {case_num}, Maekawa's theorem fails: |{n_M} - {n_V}| = {abs(n_M-n_V)}, which is not 2.")
            results.append('none')
            continue

        # Check if 't' is present. If not, can't solve for it.
        if 't' not in angles_with_t:
             results.append('none')
             continue

        # Check 2: Kawasaki's Theorem
        t_index = angles_with_t.index('t')
        
        sum_odd_known_vals = [a for idx, a in enumerate(angles_with_t) if idx % 2 == 0 and a != 't']
        sum_even_known_vals = [a for idx, a in enumerate(angles_with_t) if idx % 2 != 0 and a != 't']

        sum_odd_known = sum(sum_odd_known_vals)
        sum_even_known = sum(sum_even_known_vals)
        
        t_val = 'none'

        # 't' is at an even index (a1, a3, etc.)
        if t_index % 2 == 0:
            if round(sum_even_known) != 180:
                sum_str = " + ".join(map(str, sum_even_known_vals))
                print(f"For pattern {case_num}, Kawasaki's theorem fails: the sum of alternating known angles {sum_str} = {sum_even_known}, which is not 180.")
            else:
                t_val = 180 - sum_odd_known
                # Print the equation
                known_angles_str = " + ".join(map(str, sum_odd_known_vals))
                print(f"For pattern {case_num}, the solving equation is {known_angles_str} + t = 180.")
        # 't' is at an odd index (a2, a4, etc.)
        else: 
            if round(sum_odd_known) != 180:
                sum_str = " + ".join(map(str, sum_odd_known_vals))
                print(f"For pattern {case_num}, Kawasaki's theorem fails: the sum of alternating known angles {sum_str} = {sum_odd_known}, which is not 180.")
            else:
                t_val = 180 - sum_even_known
                # Print the equation
                known_angles_str = " + ".join(map(str, sum_even_known_vals))
                print(f"For pattern {case_num}, the solving equation is {known_angles_str} + t = 180.")
        
        if t_val != 'none' and t_val <= 0:
            results.append('none')
        elif t_val != 'none':
            results.append(int(t_val) if t_val == int(t_val) else t_val)
        else:
            results.append('none')
            
    # Print the final list of results
    print("\nFinal Result List:")
    print(f"[{','.join(map(str, results))}]")

solve_flat_foldability()