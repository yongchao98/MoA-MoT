import math

def solve_flat_foldability(pattern_data):
    """
    Analyzes a single-vertex crease pattern to find the angle 't' for flat-foldability.
    It checks both Maekawa's and Kawasaki's theorems.
    Returns a tuple containing the solution for t (or "none") and the equation used.
    """
    # Separate the list into angles and fold types
    angles_str = pattern_data[0::2]
    folds = pattern_data[1::2]

    # 1. Check Maekawa's Theorem: |#M - #V| = 2
    if abs(folds.count('M') - folds.count('V')) != 2:
        return "none", None

    # Check if 't' is in the pattern. If not, can't solve for it.
    try:
        t_index = angles_str.index('t')
    except ValueError:
        return "none", None

    # 2. Check Kawasaki's Theorem: Sum of alternating angles is 180 degrees.
    # Group angles into alternating sets (odd-positioned like a1, a3... and even-positioned like a2, a4...)
    # The list is 0-indexed, so odd-positioned angles are at even indices (0, 2, ...)
    # and even-positioned angles are at odd indices (1, 3, ...).
    
    odd_group_known_angles = [a for i, a in enumerate(angles_str) if i % 2 == 0 and a != 't']
    even_group_known_angles = [a for i, a in enumerate(angles_str) if i % 2 != 0 and a != 't']
    
    t_value = "none"
    equation = None
    
    # Case 1: 't' is in the odd-positioned group (e.g., a1, a3, ...)
    if t_index % 2 == 0:
        # The other group (even-positioned angles) must sum to 180
        if math.isclose(sum(even_group_known_angles), 180.0):
            t_value = 180.0 - sum(odd_group_known_angles)
            # Create the equation string as requested
            equation_parts = [str(a) for a in odd_group_known_angles] + ['t']
            equation = " + ".join(equation_parts) + " = 180"

    # Case 2: 't' is in the even-positioned group (e.g., a2, a4, ...)
    else:
        # The other group (odd-positioned angles) must sum to 180
        if math.isclose(sum(odd_group_known_angles), 180.0):
            t_value = 180.0 - sum(even_group_known_angles)
            # Create the equation string as requested
            equation_parts = [str(a) for a in even_group_known_angles] + ['t']
            equation = " + ".join(equation_parts) + " = 180"

    # If a solution was found, validate it (angle must be positive) and format it.
    if t_value != "none" and t_value > 0:
        # Format to integer if it's a whole number
        if math.isclose(t_value, round(t_value)):
            return int(round(t_value)), equation
        else:
            return t_value, equation
            
    # If any check fails, return "none"
    return "none", None


# Define the four crease patterns from the problem
patterns = [
    [100, 'M', 62, 'V', 22, 'M', 't', 'V', 33, 'M', 90, 'V'],
    [90, 'M', 120, 'M', 60, 'M', 90, 'M'],
    [60, 'V', 60, 'M', 120, 'M', 't', 'M'],
    [77, 'M', 15, 'M', 50, 'V', 33, 'M', 't', 'V', 130, 'M']
]

final_results = []
# Process each pattern
for i, p in enumerate(patterns):
    result, eq_str = solve_flat_foldability(p)
    final_results.append(str(result))
    if eq_str:
        # As requested, output the equation for the solvable case
        print(f"For pattern {i+1}, the solution is found using the equation:")
        print(f"  {eq_str}")
        print(f"  Solving gives t = {result}\n")

# Print the final comma-separated list of results
print(f"The final list of values for 't' is: [{','.join(final_results)}]")
