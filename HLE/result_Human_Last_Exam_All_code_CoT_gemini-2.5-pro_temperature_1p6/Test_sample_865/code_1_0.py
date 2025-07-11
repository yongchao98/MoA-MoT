import math

def solve_crease_patterns():
    """
    Analyzes four single-vertex crease patterns to find the angle 't' that
    allows for flat-folding, based on Maekawa's and Kawasaki's theorems.
    """
    patterns = [
        {'id': 1, 'data': [100,'M',62,'V',22,'M','t','V',33,'M',90,'V']},
        {'id': 2, 'data': [90,'M',120,'M',60,'M',90,'M']},
        {'id': 3, 'data': [60,'V',60,'M',120,'M','t','M']},
        {'id': 4, 'data': [77,'M',15,'M',50,'V',33,'M','t','V',130,'M']}
    ]

    final_results = []

    for p in patterns:
        case_id = p['id']
        pattern_data = p['data']
        
        print(f"Analyzing Case {case_id}: {pattern_data}")
        
        angles_raw = [item for item in pattern_data if isinstance(item, (int, str)) and item != 'M' and item != 'V']
        folds = [item for item in pattern_data if item == 'M' or item == 'V']
        
        # 1. Check Maekawa's Theorem
        num_m = folds.count('M')
        num_v = folds.count('V')
        
        if abs(num_m - num_v) != 2:
            print(f"-> Maekawa's theorem fails: |#M - #V| = |{num_m} - {num_v}| = {abs(num_m - num_v)}, which is not 2.")
            print("-> Result: none\n" + "-"*30)
            final_results.append("none")
            continue
            
        print(f"-> Maekawa's theorem holds: |#M - #V| = |{num_m} - {num_v}| = 2.")

        # 2. Find and check angles
        t_val = None
        if 't' in angles_raw:
            known_angles_sum = sum(a for a in angles_raw if a != 't')
            t_val = 360.0 - known_angles_sum
            print(f"-> From sum of angles = 360, t is calculated to be {t_val:.0f}.")
            if t_val <= 0:
                print(f"-> Angle t={t_val} must be positive. Not a valid pattern.")
                print("-> Result: none\n" + "-"*30)
                final_results.append("none")
                continue
            angles = [t_val if a == 't' else float(a) for a in angles_raw]
        else: # No 't' in pattern
            angles = [float(a) for a in angles_raw]
            if not math.isclose(sum(angles), 360.0):
                print(f"-> Sum of angles is {sum(angles)}, not 360. Invalid pattern.")
                print("-> Result: none\n" + "-"*30)
                final_results.append("none")
                continue

        # 3. Check Kawasaki's Theorem
        odd_angles = angles[0::2]
        even_angles = angles[1::2]
        sum_odd = sum(odd_angles)
        sum_even = sum(even_angles)

        print("-> Checking Kawasaki's theorem:")
        
        odd_equation = " + ".join([f"{a:.0f}" for a in odd_angles]) + f" = {sum_odd:.0f}"
        even_equation = " + ".join([f"{a:.0f}" for a in even_angles]) + f" = {sum_even:.0f}"
        
        is_flat_foldable = True
        if not math.isclose(sum_odd, 180.0):
            print(f"  - Alternating sum 1: {odd_equation} (Should be 180)")
            is_flat_foldable = False
        else:
             print(f"  - Alternating sum 1: {odd_equation} (OK)")

        if not math.isclose(sum_even, 180.0):
            print(f"  - Alternating sum 2: {even_equation} (Should be 180)")
            is_flat_foldable = False
        else:
            print(f"  - Alternating sum 2: {even_equation} (OK)")
            
        if is_flat_foldable:
            result = int(t_val) if t_val is not None else "foldable"
            if result == "foldable": # A foldable case but no t to solve for
                 final_results.append("none")
                 print("-> The pattern is flat-foldable, but no 't' was provided to solve for.")
                 print("-> Result: none\n" + "-"*30)
            else:
                 final_results.append(str(result))
                 print(f"-> Pattern is flat-foldable for t = {result}.")
                 print(f"-> Result: {result}\n" + "-"*30)
        else:
            final_results.append("none")
            print("-> Kawasaki's theorem fails.")
            print("-> Result: none\n" + "-"*30)

    print(f"Final comma-separated list: [{','.join(final_results)}]")

solve_crease_patterns()