import math

def count_flat_foldable_assignments():
    """
    Calculates the number of valid crease assignments for a partially
    assigned pattern to be flat-foldable.
    """
    # The input from the problem description
    pattern_data = [60, 'M', 30, '?', 50, '?', 70, 'V', 150, '?']

    # Step 1: Parse the input into angles and creases
    angles = [item for i, item in enumerate(pattern_data) if i % 2 == 0]
    creases = [item for i, item in enumerate(pattern_data) if i % 2 != 0]
    
    # --- Check necessary conditions for flat-foldability ---

    # Condition 1: The number of creases around the vertex must be even.
    num_creases = len(creases)
    if num_creases % 2 != 0:
        print(0)
        return

    # Condition 2: The sum of all angles must be 360 degrees.
    if sum(angles) != 360:
        print(0)
        return

    # Condition 3: Kawasaki's Theorem - sum of alternate angles must be 180 degrees.
    sum_odd = sum(angles[0::2])
    sum_even = sum(angles[1::2])
    if sum_odd != 180 or sum_even != 180:
        print(0)
        return

    # --- If geometry is valid, count assignments using Maekawa's Theorem ---
    
    # Maekawa's Theorem: The number of Mountain folds and Valley folds must differ by 2.
    # |M - V| = 2

    M_known = creases.count('M')
    V_known = creases.count('V')
    Q_count = creases.count('?')
    
    valid_assignments = 0
    
    # Iterate through all ways to assign '?' as 'M' (from 0 to Q_count)
    for m in range(Q_count + 1):
        v = Q_count - m  # The rest are assigned as 'V'
        
        M_total = M_known + m
        V_total = V_known + v
        
        if abs(M_total - V_total) == 2:
            # Number of ways to choose 'm' M-folds from 'Q_count' options
            # is the binomial coefficient "Q_count choose m".
            num_ways = math.comb(Q_count, m)
            valid_assignments += num_ways
            
    print(valid_assignments)

count_flat_foldable_assignments()