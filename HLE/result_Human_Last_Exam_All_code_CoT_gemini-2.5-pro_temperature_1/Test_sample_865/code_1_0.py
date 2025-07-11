import math

def solve_crease_pattern(pattern, case_num):
    """
    Analyzes a single-vertex crease pattern to find the angle 't' for flat-foldability.
    Prints the analysis and returns the result.
    """
    print(f"Analyzing Case {case_num}: {pattern}")

    # 1. Parse the input pattern into angles and crease types
    angles_str = [str(p) for p in pattern[0::2]]
    creases = [c for c in pattern[1::2]]

    # 2. Check Maekawa's Theorem: |#M - #V| = 2
    num_m = creases.count('M')
    num_v = creases.count('V')
    if abs(num_m - num_v) != 2:
        print(f"  - Fails Maekawa's Theorem: The difference between mountain folds ({num_m}) and valley folds ({num_v}) is {abs(num_m - num_v)}, not 2.")
        print(f"  - Result: none\n")
        return 'none'
    else:
        print(f"  - Passes Maekawa's Theorem: |#M - #V| = |{num_m} - {num_v}| = 2.")

    # Check if 't' exists in the pattern
    if 't' not in angles_str:
        # This case is for patterns that pass Maekawa's but don't have a 't'.
        # Since the problem asks to find 't', the answer is 'none' if 't' isn't present.
        print(f"  - Pattern has no unknown angle 't' to solve for.")
        print(f"  - Result: none\n")
        return 'none'
        
    t_index = angles_str.index('t')
    other_angles_values = [float(a) for i, a in enumerate(angles_str) if i != t_index]

    # 3. Use the Full Circle Condition to find 't': Sum of angles = 360
    sum_other_angles = sum(other_angles_values)
    t_sum = 360.0 - sum_other_angles
    
    sum_eq_parts = [a for a in angles_str if a != 't']
    print(f"  - Equation 1 (Sum of angles): {' + '.join(sum_eq_parts)} + t = 360")
    print(f"    {sum_other_angles} + t = 360  =>  t = {t_sum}")

    # 4. Use Kawasaki's Theorem to find 't': Sum of alternating angles are equal
    odd_angles_str = []
    even_angles_str = []
    odd_angles_sum = 0.0
    even_angles_sum = 0.0
    
    # Separate angles into odd and even positions
    for i, angle_str in enumerate(angles_str):
        if i % 2 == 0:  # Odd-numbered angles (1st, 3rd, etc.)
            odd_angles_str.append(angle_str)
            if angle_str != 't':
                odd_angles_sum += float(angle_str)
        else:  # Even-numbered angles (2nd, 4th, etc.)
            even_angles_str.append(angle_str)
            if angle_str != 't':
                even_angles_sum += float(angle_str)

    # Form the equation string
    kawasaki_eq_str = f"{' + '.join(odd_angles_str)} = {' + '.join(even_angles_str)}"
    print(f"  - Equation 2 (Kawasaki's Theorem): {kawasaki_eq_str}")
    
    # Solve for t
    t_is_in_odd_group = (t_index % 2 == 0)
    if t_is_in_odd_group:
        t_kawasaki = even_angles_sum - odd_angles_sum
        print(f"    t + {odd_angles_sum} = {even_angles_sum}  =>  t = {t_kawasaki}")
    else:
        t_kawasaki = odd_angles_sum - even_angles_sum
        print(f"    {odd_angles_sum} = t + {even_angles_sum}  =>  t = {t_kawasaki}")
        
    # 5. Compare results and determine the final answer
    if math.isclose(t_sum, t_kawasaki):
        solution = t_sum
        # Check if angle is valid (0 < t < 180)
        if 0 < solution < 180:
             final_t = int(solution) if math.isclose(solution, round(solution)) else solution
             print(f"  - The two equations give a consistent and valid angle.")
             print(f"  - Result: {final_t}\n")
             return final_t
        else:
            print(f"  - The angle t={solution} is not valid (must be between 0 and 180 degrees).")
            print(f"  - Result: none\n")
            return 'none'
    else:
        print(f"  - The equations give contradictory values for t ({t_sum} vs {t_kawasaki}).")
        print(f"  - Result: none\n")
        return 'none'

def main():
    # Define the four crease patterns
    patterns = [
        [100,'M',62,'V',22,'M','t','V',33,'M',90,'V'],
        [90,'M',120,'M',60,'M',90,'M'],
        [60,'V',60,'M',120,'M','t','M'],
        [77,'M',15,'M',50,'V',33,'M','t','V',130,'M']
    ]

    results = []
    for i, p in enumerate(patterns):
        result = solve_crease_pattern(p, i + 1)
        results.append(result)

    # Print the final list of results
    print("Final list of values for t:")
    print(results)

if __name__ == "__main__":
    main()