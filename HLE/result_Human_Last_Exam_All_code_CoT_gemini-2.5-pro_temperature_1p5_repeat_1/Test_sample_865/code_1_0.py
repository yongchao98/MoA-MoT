import math

def solve_flat_foldable_angle(pattern_data):
    """
    Calculates the unknown angle 't' that makes a single-vertex crease pattern flat-foldable.
    
    A pattern is flat-foldable if:
    1. Maekawa's Theorem: |Number of Mountain folds - Number of Valley folds| = 2.
    2. Kawasaki's Theorem: The sum of alternating angles are equal.
    3. The sum of all angles is 360 degrees.
    
    This function returns the numeric value for 't' if a unique, positive solution exists,
    otherwise, it returns 'none'.
    """
    
    # Step 1: Parse the input data into separate lists for angles and fold types.
    angles = []
    folds = []
    for i, item in enumerate(pattern_data):
        if i % 2 == 0:
            angles.append(item)
        else:
            folds.append(item)
            
    # Step 2: Check Maekawa's Theorem.
    mountain_folds = folds.count('M')
    valley_folds = folds.count('V')
    if abs(mountain_folds - valley_folds) != 2:
        return 'none'
        
    # If there is no unknown angle 't', we just check if the pattern is valid.
    # For the given problems, these cases are already ruled out by Maekawa's theorem.
    if 't' not in angles:
        return 'none'

    # Step 3: Solve for 't' using the angle sum and Kawasaki's theorem.
    t_index = -1
    known_angles_sum = 0
    for i, angle in enumerate(angles):
        if angle == 't':
            t_index = i
        else:
            known_angles_sum += angle
            
    # Equation 1: The sum of all angles must be 360.
    # sum of known angles + t = 360  =>  t = 360 - sum of known angles
    t_from_sum = 360 - known_angles_sum
    
    # Equation 2: The sum of alternating angles must be equal (Kawasaki's Theorem).
    # sum of odd-positioned angles = sum of even-positioned angles
    odd_indices_sum = 0
    even_indices_sum = 0
    for i, angle in enumerate(angles):
        if angle == 't':
            continue
        if i % 2 == 0:  # Angles at positions 1, 3, 5... (0-indexed)
            odd_indices_sum += angle
        else:           # Angles at positions 2, 4, 6... (1-indexed)
            even_indices_sum += angle
            
    t_from_kawasaki = 0
    if t_index % 2 == 0:  # If 't' is at an odd position (e.g., a1, a3...)
        t_from_kawasaki = even_indices_sum - odd_indices_sum
    else:                 # If 't' is at an even position (e.g., a2, a4...)
        t_from_kawasaki = odd_indices_sum - even_indices_sum

    # Step 4: A solution exists only if both methods yield the same positive angle.
    if math.isclose(t_from_sum, t_from_kawasaki) and t_from_sum > 0:
        return int(round(t_from_sum))
    else:
        return 'none'

def solve_all_patterns():
    """
    Solves for the unknown angle 't' for all four specified crease patterns.
    """
    patterns = [
        [100, 'M', 62, 'V', 22, 'M', 't', 'V', 33, 'M', 90, 'V'],
        [90, 'M', 120, 'M', 60, 'M', 90, 'M'],
        [60, 'V', 60, 'M', 120, 'M', 't', 'M'],
        [77, 'M', 15, 'M', 50, 'V', 33, 'M', 't', 'V', 130, 'M']
    ]
    
    results = [solve_flat_foldable_angle(p) for p in patterns]
    
    # Format the output as a comma-separated list inside brackets.
    output_string = f"[{','.join(map(str, results))}]"
    print(output_string)

# Execute the solver.
solve_all_patterns()