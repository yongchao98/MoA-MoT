import itertools

def solve_flat_fold():
    """
    Calculates the number of valid flat-foldable assignments for a given
    partially assigned crease pattern at a single vertex.
    """
    # The input crease pattern and angles
    pattern_data = [60, 'M', 30, '?', 50, '?', 70, 'V', 150, '?']

    # Step 1: Parse the input into angles and creases
    angles = [float(x) for x in pattern_data[0::2]]
    creases = [str(x) for x in pattern_data[1::2]]
    num_creases = len(creases)

    # Step 2: Check the necessary conditions for flat-foldability.

    # Condition A: The total number of creases must be even.
    # This is a consequence of the Maekawa-Justin Theorem (|N_M - N_V| = 2).
    if num_creases % 2 != 0:
        # The number of creases is odd, so it's impossible to satisfy Maekawa's theorem.
        print(0)
        return

    # Condition B: The sum of all angles must be 360 degrees.
    # We use a small tolerance for floating-point comparisons.
    if abs(sum(angles) - 360.0) > 1e-9:
        print(0)
        return

    # Condition C: Kawasaki's Theorem must hold for the angles.
    # The sum of alternating angles must be 180 degrees.
    alt_sum1 = sum(angles[0::2])
    alt_sum2 = sum(angles[1::2])
    if not (abs(alt_sum1 - 180.0) < 1e-9 and abs(alt_sum2 - 180.0) < 1e-9):
        print(0)
        return

    # Step 3: Iterate through all possible assignments for unassigned creases ('?').
    q_indices = [i for i, c in enumerate(creases) if c == '?']
    num_q = len(q_indices)
    
    valid_assignments_count = 0
    
    # Generate all combinations of 'M' (Mountain) and 'V' (Valley) for the '?' marks.
    possible_assignments = itertools.product(['M', 'V'], repeat=num_q)
    
    for assignment in possible_assignments:
        # Create a complete list of creases for the current assignment
        temp_creases = list(creases)
        for i, val in enumerate(assignment):
            temp_creases[q_indices[i]] = val
        
        # Step 4: Check Maekawa-Justin Theorem for the complete assignment.
        # The number of mountain folds must differ from valley folds by exactly 2.
        num_m = temp_creases.count('M')
        num_v = temp_creases.count('V')
        
        if abs(num_m - num_v) == 2:
            valid_assignments_count += 1
            
    # Step 5: Print the total number of valid assignments found.
    print(valid_assignments_count)

# Execute the function to find the solution.
solve_flat_fold()