import math

def solve_flat_foldable_angle(patterns):
    """
    Analyzes a list of single-vertex crease patterns to find the angle 't'
    that makes them flat-foldable.

    For a pattern to be flat-foldable, two conditions must be met:
    1. Maekawa's Theorem: |Number of Mountain Folds - Number of Valley Folds| = 2
    2. Kawasaki's Theorem: The sum of alternating angles is 180 degrees.
    """
    final_results = []
    
    for i, pattern_list in enumerate(patterns):
        print(f"--- Analyzing Pattern {i+1}: {pattern_list} ---")
        
        # 1. Parse the input
        # Note: The provided format is [angle, fold, angle, fold...].
        # We extract angles and folds into separate lists.
        angles_str = [str(pattern_list[j]) for j in range(0, len(pattern_list), 2)]
        folds = [pattern_list[j] for j in range(1, len(pattern_list), 2)]
        
        # 2. Check Maekawa's Theorem
        num_m = folds.count('M')
        num_v = folds.count('V')
        
        if abs(num_m - num_v) != 2:
            print(f"Maekawa's condition is not met.")
            print(f"Number of Mountain folds (M): {num_m}")
            print(f"Number of Valley folds (V): {num_v}")
            print(f"|{num_m} - {num_v}| = {abs(num_m - num_v)}, which is not 2.")
            print("Result: none\n")
            final_results.append('none')
            continue

        print(f"Maekawa's condition is met: |{num_m} - {num_v}| = 2.")

        # 3. Check Kawasaki's Theorem
        # No 't' present, special case
        if 't' not in angles_str:
            print("No variable 't' found in the pattern.")
            angles_num = [float(a) for a in angles_str]
            odd_angles = angles_num[0::2]
            even_angles = angles_num[1::2]
            sum_odd = sum(odd_angles)
            sum_even = sum(even_angles)
            print(f"Odd angles sum: {' + '.join(map(str, odd_angles))} = {sum_odd}")
            print(f"Even angles sum: {' + '.join(map(str, even_angles))} = {sum_even}")
            if not (math.isclose(sum_odd, 180) and math.isclose(sum_even, 180)):
                 print("Angle sums do not meet Kawasaki's theorem (both alternating sums must be 180).")
                 print("Result: none\n")
                 final_results.append('none')
            else: # This case is unlikely based on the prompt but included for completeness
                 print("The pattern is already flat-foldable.")
                 print("Result: none (no 't' to solve for)\n")
                 final_results.append('none')
            continue

        # 't' is present
        t_is_in_odd = False
        odd_angles_known = []
        even_angles_known = []
        
        for j, angle in enumerate(angles_str):
            if angle == 't':
                # j is 0-indexed, position is (j+1)
                if (j + 1) % 2 != 0: # Odd position
                    t_is_in_odd = True
            else:
                if (j + 1) % 2 != 0: # Odd position
                    odd_angles_known.append(float(angle))
                else: # Even position
                    even_angles_known.append(float(angle))
        
        # Check the sum of the list that does NOT contain 't'
        if t_is_in_odd:
            sum_even_known = sum(even_angles_known)
            print(f"Even angles (without t): {even_angles_known}")
            print(f"Sum of even angles: {' + '.join(map(str, even_angles_known))} = {sum_even_known}")
            if not math.isclose(sum_even_known, 180):
                print("This sum is not 180, so Kawasaki's theorem cannot be satisfied.")
                print("Result: none\n")
                final_results.append('none')
                continue
            
            # Solve for t using the odd angles
            sum_odd_known = sum(odd_angles_known)
            t_val = 180 - sum_odd_known
            print(f"Odd angles (with t): {odd_angles_known} and t")
            print(f"Equation from odd angles: {' + '.join(map(str, odd_angles_known))} + t = 180")
            print(f"Solving for t: t = 180 - {sum_odd_known} = {t_val}")
            print(f"Result: {t_val}\n")
            final_results.append(str(t_val))
        else: # t is in even angles
            sum_odd_known = sum(odd_angles_known)
            print(f"Odd angles (without t): {odd_angles_known}")
            print(f"Sum of odd angles: {' + '.join(map(str, odd_angles_known))} = {sum_odd_known}")
            if not math.isclose(sum_odd_known, 180):
                print("This sum is not 180, so Kawasaki's theorem cannot be satisfied.")
                print("Result: none\n")
                final_results.append('none')
                continue

            # Solve for t using the even angles
            sum_even_known = sum(even_angles_known)
            t_val = 180 - sum_even_known
            print(f"Even angles (with t): {even_angles_known} and t")
            print(f"Equation from even angles: {' + '.join(map(str, even_angles_known))} + t = 180")
            print(f"Solving for t: t = 180 - {sum_even_known} = {t_val}")
            print(f"Result: {t_val}\n")
            final_results.append(str(t_val))

    return final_results

# The four single-vertex crease patterns from the problem
patterns = [
    [100, 'M', 62, 'V', 22, 'M', 't', 'V', 33, 'M', 90, 'V'],
    [90, 'M', 120, 'M', 60, 'M', 90, 'M'],
    [60, 'V', 60, 'M', 120, 'M', 't', 'M'],
    [77, 'M', 15, 'M', 50, 'V', 33, 'M', 't', 'V', 130, 'M']
]

results = solve_flat_foldable_angle(patterns)

# Print the final summary list
print("--- Final Summary ---")
# Format the result list as requested, e.g., [72,33,none,5]
# Convert numbers to int if they are whole, otherwise float.
formatted_results = []
for res in results:
    if res == 'none':
        formatted_results.append('none')
    else:
        val = float(res)
        if val.is_integer():
            formatted_results.append(str(int(val)))
        else:
            formatted_results.append(str(val))

print(f"[{','.join(formatted_results)}]")
