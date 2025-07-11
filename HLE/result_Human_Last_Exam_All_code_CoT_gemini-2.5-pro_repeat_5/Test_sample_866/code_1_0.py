import math

def solve_flat_folding():
    """
    Calculates the number of valid flat-foldable assignments for a given
    partially assigned crease pattern at a single vertex.
    """
    # The input crease pattern
    pattern_data = [60, 'M', 30, '?', 50, '?', 70, 'V', 150, '?']

    # 1. Parse the input into angles and creases
    angles = [item for item in pattern_data if isinstance(item, (int, float))]
    creases = [item for item in pattern_data if isinstance(item, str)]
    n = len(angles)

    print(f"Analyzing crease pattern with {n} angles.")
    print(f"Angles: {angles}")
    print(f"Creases: {creases}")
    print("-" * 30)

    # 2. Check Kawasaki's Theorem
    # The sum of alternating angles must be 180 degrees.
    print("Step 1: Checking Kawasaki's Theorem...")

    # Sum of odd-indexed angles (α₁, α₃, ...)
    odd_indices_angles = [angles[i] for i in range(0, n, 2)]
    sum_odd = sum(odd_indices_angles)
    odd_sum_eq = " + ".join(map(str, odd_indices_angles))
    print(f"Sum of alternating angles (1st set): {odd_sum_eq} = {sum_odd}")

    # Sum of even-indexed angles (α₂, α₄, ...)
    even_indices_angles = [angles[i] for i in range(1, n, 2)]
    sum_even = sum(even_indices_angles)
    even_sum_eq = " + ".join(map(str, even_indices_angles))
    print(f"Sum of alternating angles (2nd set): {even_sum_eq} = {sum_even}")

    # For flat-foldability, both sums must be 180.
    # We use math.isclose for robust floating-point comparison.
    if not (math.isclose(sum_odd, 180) and math.isclose(sum_even, 180)):
        print("\nResult: Kawasaki's Theorem is NOT satisfied.")
        print("The alternating angle sums must both be 180 for the pattern to be flat-foldable.")
        print("Since the angles do not meet this condition, no assignment of creases can make it valid.")
        final_answer = 0
    else:
        # This part of the code will not be reached for the given input.
        print("\nResult: Kawasaki's Theorem is satisfied.")
        print("-" * 30)
        print("Step 2: Checking Maekawa's Theorem...")

        # Maekawa's theorem (|#M - #V| = 2) implies n must be even.
        if n % 2 != 0:
            print(f"The number of creases is {n}, which is odd.")
            print("Maekawa's Theorem cannot be satisfied for an odd number of creases.")
            final_answer = 0
        else:
            # Count known and unassigned creases
            num_m_known = creases.count('M')
            num_v_known = creases.count('V')
            num_q = creases.count('?')
            total_assignments = 0

            # Case 1: #M - #V = 2
            m_target = (n + 2) // 2
            v_target = n - m_target
            m_needed = m_target - num_m_known
            v_needed = v_target - num_v_known
            
            print(f"\nAnalyzing case #M - #V = 2 (Total M={m_target}, V={v_target})")
            if m_needed >= 0 and v_needed >= 0 and m_needed + v_needed == num_q:
                assignments1 = math.comb(num_q, m_needed)
                print(f"  Requires assigning {m_needed} 'M' and {v_needed} 'V' to the {num_q} unassigned creases.")
                print(f"  Number of ways: C({num_q}, {m_needed}) = {assignments1}")
                total_assignments += assignments1
            else:
                assignments1 = 0
                print("  This case is not possible with the given known creases.")

            # Case 2: #V - #M = 2
            m_target = (n - 2) // 2
            v_target = n - m_target
            m_needed = m_target - num_m_known
            v_needed = v_target - num_v_known

            print(f"\nAnalyzing case #V - #M = 2 (Total M={m_target}, V={v_target})")
            if m_needed >= 0 and v_needed >= 0 and m_needed + v_needed == num_q:
                assignments2 = math.comb(num_q, m_needed)
                print(f"  Requires assigning {m_needed} 'M' and {v_needed} 'V' to the {num_q} unassigned creases.")
                print(f"  Number of ways: C({num_q}, {m_needed}) = {assignments2}")
                total_assignments += assignments2
            else:
                assignments2 = 0
                print("  This case is not possible with the given known creases.")

            final_answer = total_assignments

    print("-" * 30)
    print(f"Total number of different flat-foldable assignments: {final_answer}")


# Execute the function to solve the problem
solve_flat_folding()