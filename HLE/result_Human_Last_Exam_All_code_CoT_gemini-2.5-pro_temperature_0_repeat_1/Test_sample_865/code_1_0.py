import math

def solve_crease_pattern(pattern_data, pattern_index):
    """
    Analyzes a single-vertex crease pattern to find the angle 't' for flat-foldability.
    It checks both Maekawa's and Kawasaki's theorems.
    """
    angles = []
    folds = []
    t_info = {'present': False, 'index': -1}

    # Step 1: Parse the input pattern into angles and folds
    for i in range(0, len(pattern_data), 2):
        angle = pattern_data[i]
        if angle == 't':
            angles.append(None)
            t_info['present'] = True
            t_info['index'] = len(angles) - 1
        else:
            angles.append(float(angle))
        
        if i + 1 < len(pattern_data):
            folds.append(pattern_data[i+1])

    # If 't' is not in the pattern, we cannot provide a value for it.
    if not t_info['present']:
        return 'none'

    # Step 2: Check Maekawa's Theorem: |#M - #V| = 2
    num_m = folds.count('M')
    num_v = folds.count('V')
    if abs(num_m - num_v) != 2:
        return 'none'

    # Step 3: Apply Kawasaki's Theorem to solve for 't'
    sum_odd_known = 0
    sum_even_known = 0
    odd_angles_str = []
    even_angles_str = []

    # A_1, A_3, ... are odd-indexed angles (Python index 0, 2, ...)
    # A_2, A_4, ... are even-indexed angles (Python index 1, 3, ...)
    t_is_in_odd_group = (t_info['index'] % 2 == 0)

    for i, angle in enumerate(angles):
        is_odd_group_angle = (i % 2 == 0)
        
        # Prepare strings for printing the equation later
        angle_str = 't' if angle is None else str(int(angle))
        if is_odd_group_angle:
            odd_angles_str.append(angle_str)
        else:
            even_angles_str.append(angle_str)

        # Sum the known angles
        if angle is not None:
            if is_odd_group_angle:
                sum_odd_known += angle
            else:
                sum_even_known += angle

    # Solve for t based on which group it's in
    t = 0
    if t_is_in_odd_group:
        # sum_odd_known + t = sum_even_known
        t = sum_even_known - sum_odd_known
    else:
        # sum_odd_known = sum_even_known + t
        t = sum_odd_known - sum_even_known

    # Step 4: Validate the solution for 't'
    # The angle 't' must be positive.
    if t <= 0:
        return 'none'

    # The sum of all angles must be 360 degrees.
    total_sum = sum(a for a in angles if a is not None) + t
    if not math.isclose(total_sum, 360.0):
        return 'none'
        
    # Step 5: If all conditions are met, print the equation and return the result.
    print(f"Solution found for pattern {pattern_index}:")
    if t_is_in_odd_group:
        final_odd_str = [s.replace('t', str(int(t))) for s in odd_angles_str]
        equation_str = f"{' + '.join(final_odd_str)} = {' + '.join(even_angles_str)}"
    else:
        final_even_str = [s.replace('t', str(int(t))) for s in even_angles_str]
        equation_str = f"{' + '.join(odd_angles_str)} = {' + '.join(final_even_str)}"
    
    print(equation_str)
    return int(t)

def main():
    """
    Main function to process the four crease patterns and print the final result.
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
        results.append(result)

    print(f"\nFinal list of values for t: {results}")

if __name__ == "__main__":
    main()