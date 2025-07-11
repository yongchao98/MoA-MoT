def solve_single_vertex_fold(pattern_data, case_number):
    """
    Analyzes a single-vertex crease pattern to find the angle 't'
    that allows for flat-folding.

    Args:
        pattern_data (list): A list of angles and fold types ('M', 'V', 't').
        case_number (int): The identifier for the pattern case.

    Returns:
        tuple: A tuple containing the result ('none' or the angle t)
               and a descriptive equation string if a solution is found.
    """
    angles = []
    folds = []

    # 1. Parse the input list into separate lists for angles and folds
    for item in pattern_data:
        if isinstance(item, str) and item in 'MV':
            folds.append(item)
        else:
            angles.append(item)

    # 2. Check Maekawa's Theorem: |M| - |V| = 2
    num_m = folds.count('M')
    num_v = folds.count('V')
    if abs(num_m - num_v) != 2:
        return 'none', None

    # 3. Check Kawasaki's Theorem
    t_index = -1
    if 't' in angles:
        t_index = angles.index('t')
    else:
        # If no 't' is present, the pattern is fixed.
        # We check if it's already flat-foldable, but since we are asked
        # to find a value for 't', the answer is 'none' for this case.
        return 'none', None

    # Kawasaki's Condition A: Sum of all angles is 360 degrees.
    # Solve for t.
    sum_known_angles = sum(a for a in angles if a != 't')
    t_val = 360.0 - sum_known_angles

    # An angle must be positive to be geometrically valid.
    if t_val <= 0:
        return 'none', None

    # Temporarily place the calculated t_val into the angles list for the next check.
    angles[t_index] = t_val

    # Kawasaki's Condition B: The sum of alternating angles must be equal.
    sum_odd_indexed_angles = sum(angles[i] for i in range(0, len(angles), 2))
    sum_even_indexed_angles = sum(angles[i] for i in range(1, len(angles), 2))

    # Use a small tolerance for floating point comparisons
    if abs(sum_odd_indexed_angles - sum_even_indexed_angles) > 1e-9:
        return 'none', None

    # If all checks pass, we have found a valid angle t.
    # Convert to int if it's a whole number.
    final_t = int(t_val) if t_val == int(t_val) else t_val
    angles[t_index] = final_t
    
    equation_str = f"For case {case_number}, the equation is: {' + '.join(map(str, angles))} = 360"
    return final_t, equation_str


if __name__ == "__main__":
    patterns = [
        {'id': 1, 'data': [100,'M',62,'V',22,'M','t','V',33,'M',90,'V']},
        {'id': 2, 'data': [90,'M',120,'M',60,'M',90,'M']},
        {'id': 3, 'data': [60,'V',60,'M',120,'M','t','M']},
        {'id': 4, 'data': [77,'M',15,'M',50,'V',33,'M','t','V',130,'M']}
    ]

    final_results = []
    for pattern in patterns:
        result, equation = solve_single_vertex_fold(pattern['data'], pattern['id'])
        if equation:
            print(equation)
        final_results.append(result)

    print(f"[{','.join(map(str, final_results))}]")
