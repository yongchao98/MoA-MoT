import math

def count_flat_foldable_assignments():
    """
    Calculates the number of valid crease assignments for a single-vertex
    origami pattern to be flat-foldable.
    """
    pattern = [60, 'M', 30, '?', 50, '?', 70, 'V', 150, '?']

    # Step 1: Parse the input into angles and creases
    angles = [item for i, item in enumerate(pattern) if i % 2 == 0]
    creases = [item for i, item in enumerate(pattern) if i % 2 != 0]

    n = len(creases)

    # Step 2: Check Maekawa's Theorem prerequisite: number of creases must be even.
    # If n is odd, |#M - #V| will always be an odd number and can never be 2.
    if n % 2 != 0:
        print("For a single vertex to be flat-foldable, the number of creases (n) must be even.")
        print("This is a consequence of Maekawa's Theorem, which states that the number of mountain folds and valley folds must differ by two.")
        print(f"In this pattern, the number of creases is n = {n}.")
        print("Checking the even-crease condition:")
        print(f"Final Equation: {n} % 2 == 0")
        print(f"Calculation: {n} % 2 = {n % 2}")
        print(f"Since the result is not 0, the number of creases is odd, and the condition is not met.")
        print("Therefore, no assignment of creases can make the pattern flat-foldable.")
        print(0)
        return

    # Step 3: Check Kawasaki's Theorem (only if n is even)
    sum_odd = sum(angles[i] for i in range(0, n, 2))
    sum_even = sum(angles[i] for i in range(1, n, 2))

    if not math.isclose(sum_odd, 180) or not math.isclose(sum_even, 180):
        print("For a single vertex to be flat-foldable, it must satisfy Kawasaki's Theorem.")
        print("This means the sum of alternating angles must be 180 degrees.")
        print("Checking Kawasaki's Theorem:")
        print(f"Sum of odd-indexed angles: {sum_odd}")
        print(f"Sum of even-indexed angles: {sum_even}")
        print("Since one or both sums are not 180, the condition is not met.")
        print("Therefore, there are 0 possible flat-foldable assignments.")
        print(0)
        return

    # Step 4: Count valid assignments using Maekawa's Theorem (if both prior checks pass)
    num_m_known = creases.count('M')
    num_v_known = creases.count('V')
    num_q = creases.count('?')
    
    total_assignments = 0
    
    # We need to find the number of ways to choose k '?' to be 'M' (0 <= k <= num_q)
    # such that |(num_m_known + k) - (num_v_known + num_q - k)| = 2
    
    # Case a: (num_m_known - num_v_known - num_q) + 2k = 2
    # 2k = 2 - num_m_known + num_v_known + num_q
    val1 = 2 - num_m_known + num_v_known + num_q
    if val1 % 2 == 0:
        k1 = val1 // 2
        if 0 <= k1 <= num_q:
            total_assignments += math.comb(num_q, k1)

    # Case b: (num_m_known - num_v_known - num_q) + 2k = -2
    # 2k = -2 - num_m_known + num_v_known + num_q
    val2 = -2 - num_m_known + num_v_known + num_q
    if val2 % 2 == 0:
        k2 = val2 // 2
        if 0 <= k2 <= num_q:
            total_assignments += math.comb(num_q, k2)

    print(f"The total number of valid assignments is {total_assignments}.")
    print(total_assignments)

count_flat_foldable_assignments()