import math

def count_flat_foldable_assignments():
    """
    Calculates the number of valid crease assignments for a given partially
    assigned single-vertex crease pattern to be flat-foldable.
    """
    
    # Input pattern from the user
    pattern = [60, 'M', 30, '?', 50, '?', 70, 'V', 150, '?']

    # Step 1: Parse the input into angles and creases
    angles = [item for item in pattern if isinstance(item, (int, float))]
    creases = [item for item in pattern if isinstance(item, str)]
    num_creases = len(creases)

    # Step 2: Check necessary conditions for a flat-foldable vertex.

    # Condition 1: The number of creases must be even.
    # This is a consequence of both Maekawa's and Kawasaki's theorems.
    if num_creases % 2 != 0:
        print(0)
        return

    # Condition 2: The sum of angles must be 360 degrees.
    if not math.isclose(sum(angles), 360):
        print(0)
        return

    # Condition 3: Kawasaki's Theorem. Sum of alternating angles must be equal.
    # Since the total is 360, this means each sum must be 180.
    sum_alt_angles1 = sum(angles[::2])
    if not math.isclose(sum_alt_angles1, 180):
        print(0)
        return

    # Step 3: If geometric conditions are met, count valid crease assignments.
    
    # Count known M/V creases and unassigned (?) creases.
    k = creases.count('?')
    m_known = creases.count('M')
    v_known = creases.count('V')
    
    total_valid_assignments = 0

    # From Maekawa's theorem: |(m_known + m_new) - (v_known + v_new)| = 2
    # Substitute v_new = k - m_new, which gives |(m_known - v_known - k) + 2*m_new| = 2
    constant_term = m_known - v_known - k
    
    # Case A: (m_known - v_known - k) + 2*m_new = 2
    # Solving for m_new (number of new mountain folds):
    # 2*m_new = 2 - constant_term
    if (2 - constant_term) % 2 == 0:
        m_new = (2 - constant_term) // 2
        # m_new must be a possible number of assignments (between 0 and k)
        if 0 <= m_new <= k:
            # Add the number of ways to choose m_new positions from k slots
            total_valid_assignments += math.comb(k, m_new)

    # Case B: (m_known - v_known - k) + 2*m_new = -2
    # Solving for m_new:
    # 2*m_new = -2 - constant_term
    if (-2 - constant_term) % 2 == 0:
        m_new = (-2 - constant_term) // 2
        if 0 <= m_new <= k:
            # These two cases for m_new are always distinct, so no double counting.
            total_valid_assignments += math.comb(k, m_new)
            
    print(total_valid_assignments)


# Execute the function to find the answer
count_flat_foldable_assignments()