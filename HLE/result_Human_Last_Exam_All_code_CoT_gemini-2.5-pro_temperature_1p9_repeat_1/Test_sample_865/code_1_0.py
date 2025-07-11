import math

def solve_crease_pattern():
    """
    Solves for the unknown angle 't' in a series of single-vertex crease patterns
    to make them flat-foldable.
    """
    patterns = [
        [100, 'M', 62, 'V', 22, 'M', 't', 'V', 33, 'M', 90, 'V'],
        [90, 'M', 120, 'M', 60, 'M', 90, 'M'],
        [60, 'V', 60, 'M', 120, 'M', 't', 'M'],
        [77, 'M', 15, 'M', 50, 'V', 33, 'M', 't', 'V', 130, 'M']
    ]

    final_results = []

    for i, pattern_data in enumerate(patterns):
        print(f"{i+1}) {pattern_data}")

        # 1. Parse data and check for 't'
        angles_str = [pattern_data[j] for j in range(0, len(pattern_data), 2)]
        folds = [pattern_data[j] for j in range(1, len(pattern_data), 2)]

        try:
            t_index_0based = angles_str.index('t')
        except ValueError:
            print("  Error: No unknown angle 't' found in pattern.\n")
            final_results.append('none')
            continue

        # 2. Check Maekawa's Theorem: |M - V| = 2
        num_m = folds.count('M')
        num_v = folds.count('V')
        if abs(num_m - num_v) != 2:
            print(f"  Maekawa's Theorem Fails: |M-V| = |{num_m}-{num_v}| = {abs(num_m-num_v)} != 2.")
            print("  Result: none\n")
            final_results.append('none')
            continue
        
        print(f"  Maekawa's Theorem holds: |{num_m}-{num_v}| = 2.")

        # 3. Apply Kawasaki's Theorem (combined with sum of angles = 360)
        # This implies: sum of odd angles = 180 and sum of even angles = 180.
        
        angles = []
        for a in angles_str:
            angles.append(float(a) if a != 't' else 't')

        sum_odd_known = sum(angle for idx, angle in enumerate(angles) if angle != 't' and (idx + 1) % 2 != 0)
        sum_even_known = sum(angle for idx, angle in enumerate(angles) if angle != 't' and (idx + 1) % 2 == 0)

        t_pos_is_odd = (t_index_0based + 1) % 2 != 0
        solution = 'none'

        if t_pos_is_odd:
            # 't' is an odd-indexed angle. The sum of even angles must be 180.
            known_even_angles = [str(int(a)) for a in angles if isinstance(a, float) and angles.index(a) % 2 != 0]
            print(f"  Check sum of even-indexed angles: {' + '.join(known_even_angles)} = {sum_even_known}")
            if not math.isclose(sum_even_known, 180):
                print("  The sum is not 180, so the conditions cannot be satisfied.")
            else:
                t_val = 180 - sum_odd_known
                if t_val > 0:
                    solution = int(t_val) if math.isclose(t_val, round(t_val)) else t_val
                    known_odd_angles = [str(int(a)) for a in angles if isinstance(a, float) and angles.index(a) % 2 == 0]
                    print(f"  Equation for t: {' + '.join(known_odd_angles)} + t = 180 => t = 180 - {sum_odd_known} = {solution}")
                else:
                    print(f"  Solved t = {t_val}, which is not a positive angle.")
        else: # 't' is an even-indexed angle
            # The sum of odd angles must be 180.
            known_odd_angles = [str(int(a)) for a in angles if isinstance(a, float) and angles.index(a) % 2 == 0]
            print(f"  Check sum of odd-indexed angles: {' + '.join(known_odd_angles)} = {sum_odd_known}")
            if not math.isclose(sum_odd_known, 180):
                print("  The sum is not 180, so the conditions cannot be satisfied.")
            else:
                t_val = 180 - sum_even_known
                if t_val > 0:
                    solution = int(t_val) if math.isclose(t_val, round(t_val)) else t_val
                    known_even_angles = [str(int(a)) for a in angles if isinstance(a, float) and angles.index(a) % 2 != 0]
                    print(f"  Equation for t: {' + '.join(known_even_angles)} + t = 180 => t = 180 - {sum_even_known} = {solution}")
                else:
                    print(f"  Solved t = {t_val}, which is not a positive angle.")
        
        print(f"  Result: {solution}\n")
        final_results.append(solution)
        
    print("----------------------------------------")
    print("Final comma-separated list:")
    print(f"[{','.join(map(str, final_results))}]")


solve_crease_pattern()
<<<[none,none,120,none]>>>