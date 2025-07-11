import math

def calculate_period(n, m):
    """Calculates the period metric P = n^2 + m^2 + nm."""
    return n**2 + m**2 + n * m

def solve_periods():
    """
    Calculates the periods for the four given sequences of points on the hexagonal grid.
    """
    
    # Map points to their primary hexagon center's lattice coordinates (n, m)
    # where C = n*A + m*B.
    # A = C23 - C13, B = C31 - C13.
    # C13: (0,0), C23: (1,0), C31: (0,1)
    point_to_center_map = {}
    for i in range(1, 14): point_to_center_map[i] = (0, 0)
    for i in range(14, 24): point_to_center_map[i] = (1, 0)
    for i in range(24, 32): point_to_center_map[i] = (0, 1)

    # Define points on shared boundaries and their possible centers.
    shared_points_map = {
        7: [(0, 0), (1, 0)],   # Shared by H13, H23
        8: [(0, 0), (0, 1)],   # Shared by H13, H31
        22: [(0, 0), (1, 0)],  # Shared by H13, H23
        30: [(0, 0), (0, 1)],  # Shared by H13, H31
    }
    
    sequences = [
        [13, 31, 23],
        [10, 4, 23, 31],
        [5, 15, 17, 19, 21, 7],
        [4, 5, 14, 23, 18, 19, 21, 22, 31, 30, 9, 8, 13]
    ]

    final_periods = []

    for i, seq in enumerate(sequences):
        unambiguous_sum_n, unambiguous_sum_m = 0, 0
        shared_in_seq = []

        for p in seq:
            if p in shared_points_map:
                shared_in_seq.append(p)
            else:
                n, m = point_to_center_map[p]
                unambiguous_sum_n += n
                unambiguous_sum_m += m
        
        if not shared_in_seq:
            period = calculate_period(unambiguous_sum_n, unambiguous_sum_m)
            final_periods.append(period)
            continue

        min_period = float('inf')
        
        # Iterate through all assignment combinations for shared points.
        num_shared = len(shared_in_seq)
        for j in range(2**num_shared):
            temp_n, temp_m = unambiguous_sum_n, unambiguous_sum_m
            
            temp_j = j
            for k in range(num_shared):
                point = shared_in_seq[k]
                choice_index = temp_j % 2
                temp_j //= 2
                
                n_choice, m_choice = shared_points_map[point][choice_index]
                temp_n += n_choice
                temp_m += m_choice

            current_period = calculate_period(temp_n, temp_m)
            if current_period < min_period:
                min_period = current_period
        
        final_periods.append(min_period)

    print(','.join(map(str, final_periods)))

solve_periods()