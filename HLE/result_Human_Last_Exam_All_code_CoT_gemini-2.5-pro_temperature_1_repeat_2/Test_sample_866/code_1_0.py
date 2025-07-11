import math

def solve_flat_folding_assignments():
    """
    Calculates the total number of different assignments of creases
    which will make the given pattern flat-foldable.
    """
    # The partially assigned crease pattern given in the problem
    pattern_data = [60, 'M', 30, '?', 50, '?', 70, 'V', 150, '?']

    # Step 1: Parse the input data into separate lists for angles and creases.
    angles = pattern_data[0::2]
    creases = pattern_data[1::2]

    # Step 2: Check the necessary conditions for flat-foldability.
    
    # Condition 1: Maekawa's Theorem implies the number of creases must be even.
    # Let M be the number of mountain folds and V be the number of valley folds.
    # The theorem states |M - V| = 2.
    # The total number of creases N = M + V.
    # If M - V = 2, then M = V + 2, so N = (V + 2) + V = 2V + 2 (even).
    # If V - M = 2, then V = M + 2, so N = M + (M + 2) = 2M + 2 (even).
    # In both cases, N must be even.
    n_creases = len(creases)
    if n_creases % 2 != 0:
        # The number of creases is 5, which is odd. The pattern can never be
        # flat-foldable.
        print(0)
        return

    # Condition 2: Kawasaki's Theorem states the alternating sum of angles
    # must be zero for the vertex to fold flat. This is a geometric constraint.
    alternating_sum = sum(angle * ((-1)**i) for i, angle in enumerate(angles))
    if not math.isclose(alternating_sum, 0):
        # If the alternating sum is not zero, no assignment can make it fold flat.
        print(0)
        return

    # If both geometric conditions are met, proceed to count valid assignments.
    # (Note: For the given input, the code will have already exited.)
    
    M_known = creases.count('M')
    V_known = creases.count('V')
    Q_known = creases.count('?')
    
    total_valid_assignments = 0

    # We need to assign 'm' mountains and 'v' valleys to the '?'s, where m+v=Q_known.
    # Find the number of 'M' assignments (m) that satisfy Maekawa's theorem.

    # Case A: Total M - Total V = 2
    # (M_known + m) - (V_known + (Q_known - m)) = 2
    # 2m = 2 - M_known + V_known + Q_known
    numerator1 = 2 - M_known + V_known + Q_known
    if numerator1 >= 0 and numerator1 % 2 == 0:
        m1 = numerator1 // 2
        if 0 <= m1 <= Q_known:
            total_valid_assignments += math.comb(Q_known, m1)

    # Case B: Total V - Total M = 2 (or Total M - Total V = -2)
    # (M_known + m) - (V_known + (Q_known - m)) = -2
    # 2m = -2 - M_known + V_known + Q_known
    numerator2 = -2 - M_known + V_known + Q_known
    if numerator2 >= 0 and numerator2 % 2 == 0:
        m2 = numerator2 // 2
        if 0 <= m2 <= Q_known:
            total_valid_assignments += math.comb(Q_known, m2)
            
    print(total_valid_assignments)

solve_flat_folding_assignments()