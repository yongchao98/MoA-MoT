import math

def count_combinations(n, k):
    """Calculates the number of combinations 'n choose k'."""
    if k < 0 or k > n:
        return 0
    # math.comb is available in Python 3.8+ and is preferred for its accuracy and efficiency.
    try:
        return math.comb(n, k)
    except AttributeError:
        # Fallback for Python versions older than 3.8
        if k == 0 or k == n:
            return 1
        if k > n // 2:
            k = n - k
        
        res = 1
        for i in range(k):
            res = res * (n - i) // (i + 1)
        return res

def solve_crease_pattern():
    """
    Analyzes a partially assigned crease pattern to find the number of valid
    flat-foldable assignments.
    """
    pattern = [60, 'M', 30, '?', 50, '?', 70, 'V', 150, '?']

    # 1. Parse the input into separate lists for angles and creases.
    angles = [item for item in pattern if isinstance(item, (int, float))]
    creases = [item for item in pattern if isinstance(item, str)]
    n_creases = len(angles)

    print(f"Analyzing crease pattern with {n_creases} angles and {n_creases} creases.")
    print("-" * 30)

    # 2. Check if the number of creases is even.
    if n_creases % 2 != 0:
        print("Condition 1: Number of creases must be even.")
        print(f"Result: FAILED. The pattern has {n_creases} creases, which is an odd number.")
        print("\nA single vertex with an odd number of creases cannot be folded flat.")
        print("Total number of valid assignments: 0")
        print("<<<0>>>")
        return

    # 3. Check Kawasaki's Theorem for angles.
    print("Condition 2: Sum of alternate angles must be 180 degrees (Kawasaki's Theorem).")
    sum_odd_angles = sum(angles[i] for i in range(0, n_creases, 2))
    sum_even_angles = sum(angles[i] for i in range(1, n_creases, 2))
    
    # Print the equation for the sum of alternating angles
    odd_angle_strs = [str(a) for a in angles[0::2]]
    print(f"Sum of 1st, 3rd, ... angles: {' + '.join(odd_angle_strs)} = {sum_odd_angles}")
    even_angle_strs = [str(a) for a in angles[1::2]]
    print(f"Sum of 2nd, 4th, ... angles: {' + '.join(even_angle_strs)} = {sum_even_angles}")

    if not (math.isclose(sum_odd_angles, 180) and math.isclose(sum_even_angles, 180)):
        print("Result: FAILED. The sums are not both 180.")
        print("\nSince the angle condition is not met, the pattern cannot be folded flat.")
        print("Total number of valid assignments: 0")
        print("<<<0>>>")
        return
    
    print("Result: PASSED. The angle sums are correct.")
    print("-" * 30)

    # 4. Check Maekawa's Theorem for crease assignments.
    print("Condition 3: Number of Mountain (M) and Valley (V) folds must differ by 2 (|#M - #V| = 2).")
    n_M_known = creases.count('M')
    n_V_known = creases.count('V')
    n_unknown = creases.count('?')
    
    print(f"The pattern has {n_M_known} assigned M, {n_V_known} assigned V, and {n_unknown} unassigned '?' folds.")
    print(f"Let k be the number of '?' assigned as M. Then {n_unknown}-k will be assigned as V.")

    total_assignments = 0

    # Case 1: #M - #V = 2
    # Equation: (n_M_known + k) - (n_V_known + n_unknown - k) = 2
    # Solved for k: 2k = 2 - n_M_known + n_V_known + n_unknown
    print("\nSolving for k in case #M - #V = 2:")
    numerator1 = 2 - n_M_known + n_V_known + n_unknown
    print(f"  2*k = 2 - {n_M_known} + {n_V_known} + {n_unknown}")
    print(f"  2*k = {numerator1}")

    if numerator1 >= 0 and numerator1 % 2 == 0:
        k1 = numerator1 // 2
        print(f"  k = {k1}")
        if 0 <= k1 <= n_unknown:
            combinations1 = count_combinations(n_unknown, k1)
            print(f"  This requires choosing {k1} of the {n_unknown} unassigned folds to be M.")
            print(f"  Number of ways: C({n_unknown}, {k1}) = {combinations1}")
            total_assignments += combinations1
        else:
            print(f"  This value of k is invalid as it is not in the range [0, {n_unknown}].")
    else:
        print("  This equation yields no valid integer solution for k.")

    # Case 2: #V - #M = 2
    # Equation: (n_V_known + n_unknown - k) - (n_M_known + k) = 2
    # Solved for k: 2k = n_V_known - n_M_known + n_unknown - 2
    print("\nSolving for k in case #V - #M = 2:")
    numerator2 = n_V_known - n_M_known + n_unknown - 2
    print(f"  2*k = {n_V_known} - {n_M_known} + {n_unknown} - 2")
    print(f"  2*k = {numerator2}")

    if numerator2 >= 0 and numerator2 % 2 == 0:
        k2 = numerator2 // 2
        print(f"  k = {k2}")
        if 0 <= k2 <= n_unknown:
            combinations2 = count_combinations(n_unknown, k2)
            print(f"  This requires choosing {k2} of the {n_unknown} unassigned folds to be M.")
            print(f"  Number of ways: C({n_unknown}, {k2}) = {combinations2}")
            total_assignments += combinations2
        else:
            print(f"  This value of k is invalid as it is not in the range [0, {n_unknown}].")
    else:
        print("  This equation yields no valid integer solution for k.")
    
    print("-" * 30)
    print(f"The total number of valid assignments is {total_assignments}.")
    print(f"<<<{total_assignments}>>>")

if __name__ == '__main__':
    solve_crease_pattern()