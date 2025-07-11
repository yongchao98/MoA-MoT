def solve_crease_pattern():
    """
    Analyzes a partially assigned crease pattern to find the number of valid
    flat-foldable assignments.
    """
    pattern = [60, 'M', 30, '?', 50, '?', 70, 'V', 150, '?']

    # Extract angles and creases from the pattern
    angles = [item for item in pattern if isinstance(item, int)]

    print("To determine the number of valid flat-foldable assignments, we must check two conditions.")
    print("First, we check the geometric condition: Kawasaki's Theorem.\n")
    print("Kawasaki's Theorem states that the sums of alternating angles around the vertex must be equal.")

    # Separate angles into two alternating sets
    odd_indexed_angles = angles[0::2]
    even_indexed_angles = angles[1::2]

    # Calculate the sums
    sum1 = sum(odd_indexed_angles)
    sum2 = sum(even_indexed_angles)

    # Create the equation strings for printing
    sum1_eq_str = " + ".join(map(str, odd_indexed_angles))
    sum2_eq_str = " + ".join(map(str, even_indexed_angles))

    print(f"The given angles are: {angles}")
    print("Let's calculate the two alternating sums:")
    print(f"Sum 1 (1st, 3rd, 5th... angles): {sum1_eq_str} = {sum1}")
    print(f"Sum 2 (2nd, 4th, 6th... angles): {sum2_eq_str} = {sum2}\n")

    # Check if the sums are equal
    if sum1 != sum2:
        print(f"Result: The sums are not equal ({sum1} != {sum2}).")
        print("Since Kawasaki's Theorem is not satisfied, the pattern is not geometrically flat-foldable.")
        print("Therefore, no assignment of Mountain or Valley creases can make it valid.")
        final_count = 0
    else:
        # This part of the code would run if Kawasaki's theorem was satisfied.
        # We would then check Maekawa's Theorem.
        print("Result: The sums are equal. Kawasaki's Theorem is satisfied.")
        print("Now, we check the combinatorial condition: Maekawa's Theorem (|M-V|=2).\n")
        creases = [item for item in pattern if isinstance(item, str)]
        num_m_known = creases.count('M')
        num_v_known = creases.count('V')
        num_q = creases.count('?')
        
        final_count = 0
        # Iterate through possibilities for '?' (k = number of '?' assigned as 'M')
        for k in range(num_q + 1):
            total_m = num_m_known + k
            total_v = num_v_known + (num_q - k)
            if abs(total_m - total_v) == 2:
                # The number of ways to choose k items from num_q is C(num_q, k)
                import math
                combinations = math.comb(num_q, k)
                final_count += combinations

    print(f"\nThe total number of different assignments that make the pattern flat-fordable is: {final_count}")

solve_crease_pattern()
<<<0>>>