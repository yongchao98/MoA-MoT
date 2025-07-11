import math

def solve_flat_fold():
    """
    This function calculates the number of valid flat-foldable assignments
    for a given partially assigned crease pattern at a single vertex.
    """
    # The input pattern as a list of angles and crease types.
    pattern = [60, 'M', 30, '?', 50, '?', 70, 'V', 150, '?']

    # Step 1: Parse the input into separate lists for angles and creases.
    angles = [item for item in pattern if isinstance(item, int)]
    
    print("Step 1: Checking Kawasaki's Theorem")
    print("The sum of alternating angles around a vertex must be 180 degrees for a flat fold.")
    print("-" * 70)

    # Step 2: Verify Kawasaki's Theorem.
    # The problem describes a sequence of angles around a vertex.
    # For a flat fold, the sum of odd-numbered angles and the sum of
    # even-numbered angles must both be 180 degrees.
    # We use 0-based indexing for our list, so we sum elements at
    # even indices (1st, 3rd, 5th...) and odd indices (2nd, 4th, 6th...).
    
    odd_positioned_angles = angles[0::2]
    even_positioned_angles = angles[1::2]
    
    sum_odd = sum(odd_positioned_angles)
    sum_even = sum(even_positioned_angles)

    # Print the equations for the sums, showing each number as requested.
    odd_sum_eq = " + ".join(map(str, odd_positioned_angles))
    even_sum_eq = " + ".join(map(str, even_positioned_angles))
    
    print(f"Sum of 1st, 3rd, 5th angles: {odd_sum_eq} = {sum_odd}")
    print(f"Sum of 2nd, 4th angles: {even_sum_eq} = {sum_even}")
    
    print("-" * 70)

    # Step 3: Conclude based on Kawasaki's Theorem.
    # If the sums are not 180, the geometry is not flat-foldable.
    if sum_odd != 180 or sum_even != 180:
        print("Result: Kawasaki's Theorem is NOT satisfied because the alternating angle sums are not equal to 180.")
        print("A pattern that does not satisfy Kawasaki's theorem cannot be folded flat, regardless of crease assignments.")
        final_count = 0
    else:
        # This part of the logic would be executed if Kawasaki's Theorem were satisfied.
        # It uses Maekawa's Theorem to count valid assignments.
        print("Result: Kawasaki's Theorem is satisfied.")
        print("\nStep 2: Checking Maekawa's Theorem")
        print("The number of Mountain folds and Valley folds must differ by 2 ( |#M - #V| = 2 ).")
        print("-" * 70)

        creases = [item for item in pattern if isinstance(item, str)]
        num_m_known = creases.count('M')
        num_v_known = creases.count('V')
        num_q = creases.count('?')
        
        total_assignments = 0
        
        # We need to find the number of ways (m_q) to assign '?' to be M-folds.
        # Let m_q be the number of '?' assigned as M. Then (num_q - m_q) are assigned as V.
        # The total M folds will be M_total = num_m_known + m_q
        # The total V folds will be V_total = num_v_known + (num_q - m_q)

        # Case 1: M_total - V_total = 2
        # (num_m_known + m_q) - (num_v_known + num_q - m_q) = 2
        # num_m_known - num_v_known - num_q + 2*m_q = 2
        # m_q_num = (2 - num_m_known + num_v_known + num_q)
        if (2 - num_m_known + num_v_known + num_q) % 2 == 0:
            m_q1 = (2 - num_m_known + num_v_known + num_q) // 2
            if 0 <= m_q1 <= num_q:
                combinations1 = math.comb(num_q, m_q1)
                total_assignments += combinations1
        
        # Case 2: M_total - V_total = -2
        # (num_m_known + m_q) - (num_v_known + num_q - m_q) = -2
        # m_q_num = (-2 - num_m_known + num_v_known + num_q)
        if (-2 - num_m_known + num_v_known + num_q) % 2 == 0:
            m_q2 = (-2 - num_m_known + num_v_known + num_q) // 2
            if 0 <= m_q2 <= num_q and m_q1 != m_q2: # check m_q1 != m_q2 to avoid double counting
                combinations2 = math.comb(num_q, m_q2)
                total_assignments += combinations2

        final_count = total_assignments

    print("\nFinal Answer:")
    print(f"The total number of different flat-foldable assignments is: {final_count}")

solve_flat_fold()
<<<0>>>