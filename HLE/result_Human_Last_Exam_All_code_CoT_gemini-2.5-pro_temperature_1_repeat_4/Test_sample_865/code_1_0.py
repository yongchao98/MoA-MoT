import math

def solve_crease_pattern(pattern_data, case_num):
    """
    Analyzes a single-vertex crease pattern to find the angle 't' for flat-foldability.
    """
    print(f"--- Analyzing Case {case_num} ---")
    
    # 1. Parse the input data
    angles = [item for item in pattern_data if not isinstance(item, str) or item != 't']
    folds = [item for item in pattern_data if isinstance(item, str) and item in ('M', 'V')]
    
    has_t = 't' in pattern_data
    if has_t:
        t_index = pattern_data.index('t') // 2
        # The fold type for 't' is the next item in the list
        t_fold_type = pattern_data[pattern_data.index('t') + 1]
    
    # 2. Check Maekawa's Theorem: |#M - #V| = 2
    num_m = folds.count('M')
    num_v = folds.count('V')

    if abs(num_m - num_v) != 2:
        print(f"Result: Not flat-foldable.")
        print(f"Reason: Maekawa's theorem is not satisfied. The pattern has {num_m} mountain folds and {num_v} valley folds. The difference must be 2, but it is {abs(num_m - num_v)}.")
        return 'none'
        
    # If there's no 't' and Maekawa's holds, check Kawasaki's for the fixed pattern
    if not has_t:
        print(f"Result: Not flat-foldable.")
        print(f"Reason: The pattern does not contain a variable angle 't' to solve for, and it is not flat-foldable as given.")
        # We can also show why Kawasaki's fails for this specific case
        sum_odd = sum(angles[i] for i in range(0, len(angles), 2))
        sum_even = sum(angles[i] for i in range(1, len(angles), 2))
        print(f"The alternating angle sums are {sum_odd} and {sum_even}, which are not both equal to 180.")
        return 'none'

    # 3. Apply Kawasaki's Theorem: Sum of alternating angles = 180
    odd_angles_vals = []
    odd_angles_str = []
    even_angles_vals = []
    even_angles_str = []
    
    # Re-parse angles including 't' as a placeholder
    all_angles = [pattern_data[i] for i in range(0, len(pattern_data), 2)]

    for i, angle in enumerate(all_angles):
        if (i + 1) % 2 != 0: # Odd positions (1, 3, 5...)
            odd_angles_str.append(str(angle))
            if angle != 't':
                odd_angles_vals.append(angle)
        else: # Even positions (2, 4, 6...)
            even_angles_str.append(str(angle))
            if angle != 't':
                even_angles_vals.append(angle)

    # Determine which sum contains 't' and which is fully known
    sum_known, sum_unknown_parts, str_known, str_unknown = (0, [], [], [])
    if 't' in odd_angles_str:
        sum_known = sum(even_angles_vals)
        sum_unknown_parts = odd_angles_vals
        str_known = " + ".join(even_angles_str)
        str_unknown = " + ".join(odd_angles_str)
    else:
        sum_known = sum(odd_angles_vals)
        sum_unknown_parts = even_angles_vals
        str_known = " + ".join(odd_angles_str)
        str_unknown = " + ".join(even_angles_str)

    # Check if the known sum is 180
    if not math.isclose(sum_known, 180.0):
        print(f"Result: Not flat-foldable.")
        print(f"Reason: Kawasaki's theorem is not satisfied. The sum of one set of alternating angles must be 180, but it is {str_known} = {sum_known}.")
        return 'none'
    
    # Solve for t
    t_val = 180 - sum(sum_unknown_parts)
    
    # Check for valid angle
    if t_val <= 0 or t_val >= 360:
         print(f"Result: Not flat-foldable.")
         print(f"Reason: Solving for t yields an invalid angle value of {t_val}.")
         return 'none'
         
    print("Result: Solvable.")
    print("Reason: The pattern satisfies Maekawa's and Kawasaki's theorems.")
    print(f"The two alternating angle sums must be 180 degrees.")
    print(f"Sum 1: {str_known} = {sum_known}")
    print(f"Sum 2: {str_unknown} = 180")
    print(f"Solving for t: t = 180 - {' - '.join(map(str, sum_unknown_parts))}")
    print(f"Final Answer: t = {t_val}")
    return t_val


def main():
    """
    Main function to run the analysis for all four cases.
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
        if isinstance(result, (int, float)):
             # Convert to int if it's a whole number
            results.append(int(result) if result == int(result) else result)
        else:
            results.append(result)
        print()

    # Format the final list as a string
    final_list_str = f"[{','.join(map(str, results))}]"
    print("--- Summary ---")
    print("The values for t that make each crease pattern flat foldable are:")
    print(final_list_str)

if __name__ == "__main__":
    main()