import math

def solve_crease_pattern(pattern_data, case_number):
    """
    Solves for the unknown angle 't' in a single-vertex crease pattern
    to make it flat-foldable.
    """
    print(f"--- Solving Case {case_number} ---")
    print(f"Pattern: {pattern_data}")

    # Step 1: Parse angles and folds
    angles_with_t = [item for i, item in enumerate(pattern_data) if i % 2 == 0]
    folds = [item for i, item in enumerate(pattern_data) if i % 2 != 0]

    # Step 2: Check Maekawa's Theorem: |#M - #V| = 2
    num_m = folds.count('M')
    num_v = folds.count('V')
    print(f"Found {num_m} Mountain (M) folds and {num_v} Valley (V) folds.")
    if abs(num_m - num_v) != 2:
        print(f"Failed Maekawa's Theorem: |{num_m} - {num_v}| = {abs(num_m - num_v)}, which is not 2.")
        print("Result: none\n")
        return "none"
    
    print("Passed Maekawa's Theorem: |#M - #V| = 2.")

    # Check if 't' exists. If not, it can't be solved for.
    if 't' not in angles_with_t:
        print("The pattern does not contain an unknown angle 't' to solve for.")
        print("Result: none\n")
        return "none"

    # Step 3: Apply Kawasaki's Theorem
    # Replace 't' with a placeholder (None) for calculations
    angles = []
    for angle in angles_with_t:
        if angle == 't':
            angles.append(None)
        else:
            angles.append(float(angle))

    odd_angles = angles[0::2]
    even_angles = angles[1::2]
    
    # Check the group of angles that *does not* contain 't'
    if None in odd_angles:
        group_to_check = even_angles
        group_with_t = odd_angles
        group_name = "even-indexed"
    else:
        group_to_check = odd_angles
        group_with_t = even_angles
        group_name = "odd-indexed"

    sum_of_known_group = sum(group_to_check)
    known_group_str = " + ".join(map(str, [int(x) for x in group_to_check]))

    if not math.isclose(sum_of_known_group, 180.0):
        print(f"Failed Kawasaki's Theorem: The sum of {group_name} angles ({known_group_str}) is {sum_of_known_group}, not 180.")
        print("Result: none\n")
        return "none"
    
    print(f"Passed Kawasaki's Theorem check for known angles: {known_group_str} = {int(sum_of_known_group)}.")

    # Step 4: Solve for 't' in the other group
    known_angles_in_t_group = [a for a in group_with_t if a is not None]
    sum_of_knowns_in_t_group = sum(known_angles_in_t_group)
    
    t = 180.0 - sum_of_knowns_in_t_group
    
    equation_known_parts = [str(int(a)) for a in known_angles_in_t_group]
    initial_equation_str = " + ".join(equation_known_parts) + " + t = 180"
    print(f"Solving equation for 't': {initial_equation_str}")
    
    # Output the final equation with the calculated value of t
    final_eq_parts = equation_known_parts + [str(round(t, 2))]
    final_equation_str = " + ".join(final_eq_parts)
    print(f"Final Equation: {final_equation_str} = 180")
    
    # Check if the calculated angle is valid (0 < t < 180)
    if t <= 0 or t >= 180:
        print(f"Calculated value t={t} is not a valid angle.")
        print("Result: none\n")
        return "none"

    print(f"Result: {round(t, 2)}\n")
    return round(t, 2)


# --- Main execution ---

# Define the four crease patterns
patterns = [
    [100, 'M', 62, 'V', 22, 'M', 't', 'V', 33, 'M', 90, 'V'],
    [90, 'M', 120, 'M', 60, 'M', 90, 'M'],
    [60, 'V', 60, 'M', 120, 'M', 't', 'M'],
    [77, 'M', 15, 'M', 50, 'V', 33, 'M', 't', 'V', 130, 'M']
]

# Solve each pattern and store the result
results = []
for i, p in enumerate(patterns):
    result = solve_crease_pattern(p, i + 1)
    results.append(result)

# Format the final list for printing
final_output_list = [str(r) if isinstance(r, (int, float)) else "none" for r in results]

print("---")
print("Final combined answer:")
print(f"[{','.join(final_output_list)}]")
print("<<<[none,none,120,none]>>>")
