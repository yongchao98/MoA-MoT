import math

def solve_flat_foldable(pattern_list):
    """
    Analyzes a single-vertex crease pattern to find the unknown angle 't'
    that makes it flat-foldable based on Maekawa's and Kawasaki's theorems.

    Args:
        pattern_list: A list representing the crease pattern,
                      e.g., [100, 'M', 62, 'V', 't', 'M', ...].

    Returns:
        The value of the angle 't' (as an int if whole) if a solution exists, 
        otherwise the string 'none'.
    """
    angles = []
    folds = []
    t_present = False
    t_index = -1

    # Step 1: Parse the input list into angles and folds
    for i, item in enumerate(pattern_list):
        if i % 2 == 0:  # Angle value
            if item == 't':
                angles.append('t')
                t_present = True
                t_index = len(angles) - 1
            else:
                angles.append(float(item))
        else:  # Fold type
            folds.append(item)

    # Step 2: Check Maekawa's Theorem: |#M - #V| = 2
    num_m = folds.count('M')
    num_v = folds.count('V')
    if abs(num_m - num_v) != 2:
        return 'none'
    
    # The number of creases must be even for flat folding.
    if len(angles) % 2 != 0:
        return 'none'

    # Step 3 & 4: Apply Kawasaki's Theorem and solve for t
    # If 't' is not in the pattern, check if it's valid. If so, there's no 't' to solve for.
    if not t_present:
        if not math.isclose(sum(angles), 360.0) or \
           not math.isclose(sum(angles[0::2]), sum(angles[1::2])):
            return 'none' # The fixed pattern is not flat-foldable
        return 'none' # The fixed pattern is flat-foldable, but no 't' to report.

    # Calculate 't' from the sum of angles condition
    known_angles_sum = sum(angle for angle in angles if isinstance(angle, float))
    t_from_sum = 360.0 - known_angles_sum

    # Calculate 't' from the alternating angles condition
    sum_alt_even_indices = sum(angles[i] for i in range(0, len(angles), 2) if angles[i] != 't')
    sum_alt_odd_indices = sum(angles[i] for i in range(1, len(angles), 2) if angles[i] != 't')

    t_from_alternating = 0
    if t_index % 2 == 0:  # 't' is at an even index (part of the first alternating sum)
        t_from_alternating = sum_alt_odd_indices - sum_alt_even_indices
    else:  # 't' is at an odd index (part of the second alternating sum)
        t_from_alternating = sum_alt_even_indices - sum_alt_odd_indices

    # The two calculated values for 't' must be consistent (equal) and positive
    if math.isclose(t_from_sum, t_from_alternating) and t_from_sum > 0:
        # Return as an integer if it's a whole number
        if math.isclose(t_from_sum, round(t_from_sum)):
            return int(round(t_from_sum))
        return t_from_sum
    else:
        return 'none'

def main():
    """
    Main function to solve the four given crease patterns and print the result.
    """
    patterns = [
        [100, 'M', 62, 'V', 22, 'M', 't', 'V', 33, 'M', 90, 'V'],
        [90, 'M', 120, 'M', 60, 'M', 90, 'M'],
        [60, 'V', 60, 'M', 120, 'M', 't', 'M'],
        [77, 'M', 15, 'M', 50, 'V', 33, 'M', 't', 'V', 130, 'M']
    ]

    solutions = [solve_flat_foldable(p) for p in patterns]

    # Step 5: Print the consolidated results in the required list format
    print(f"[{','.join(map(str, solutions))}]")

if __name__ == "__main__":
    main()