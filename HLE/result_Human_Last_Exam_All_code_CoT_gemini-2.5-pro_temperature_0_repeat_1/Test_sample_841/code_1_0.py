def solve_unique_length_factorizations():
    """
    This function calculates the size of the set of specific quadratic integer rings
    that are Half-Factorial Domains (HFDs), meaning prime factorizations have unique lengths.
    """
    # Known lists of square-free d > 0 for which the class number h of Q(sqrt(-d)) is 1 or 2.
    h1_d = {1, 2, 3, 7, 11, 19, 43, 67, 163}
    h2_d = {5, 6, 10, 13, 15, 22, 35, 37, 51, 58, 91, 115, 123, 187, 235, 267, 403, 427}
    
    print("The problem asks for the number of rings with unique factorization lengths (HFDs).")
    print("This is equivalent to finding rings with class number 1 or 2.")
    print("The set of rings can be partitioned into two disjoint types:")
    print("Type 1: Rings of the form Z[sqrt(-d)]")
    print("Type 2: Rings of the form Z[(1+sqrt(-d))/2]\n")

    # --- Count HFDs of Type 1: Z[sqrt(-d)] ---
    print("--- Counting HFDs of Type 1: Z[sqrt(-d)] ---")
    
    # Subcase 1.1: Maximal orders (d = 1, 2 mod 4)
    # These are HFDs if h(Q(sqrt(-d))) is 1 or 2.
    d_max_h1 = sorted([d for d in h1_d if d % 4 in [1, 2]])
    d_max_h2 = sorted([d for d in h2_d if d % 4 in [1, 2]])
    count_max = len(d_max_h1) + len(d_max_h2)
    print(f"Found {count_max} HFDs that are maximal orders of this type.")
    print(f"  d values (h=1): {d_max_h1}")
    print(f"  d values (h=2): {d_max_h2}")

    # Subcase 1.2: Non-maximal orders (d = 3 mod 4)
    # h(Z[sqrt(-d)]) is 1 or 2 only for d = 3, 7, 15.
    d_non_max = [3, 7, 15]
    count_non_max = len(d_non_max)
    print(f"\nFound {count_non_max} HFDs that are non-maximal orders of this type.")
    print(f"  d values: {d_non_max}")

    count1 = count_max + count_non_max
    print(f"\nTotal HFDs of Type 1 = {count_max} + {count_non_max} = {count1}\n")

    # --- Count HFDs of Type 2: Z[(1+sqrt(-d))/2] ---
    print("--- Counting HFDs of Type 2: Z[(1+sqrt(-d))/2] ---")
    # These are maximal orders for d = 3 mod 4.
    # They are HFDs if h(Q(sqrt(-d))) is 1 or 2.
    d_type2_h1 = sorted([d for d in h1_d if d % 4 == 3])
    d_type2_h2 = sorted([d for d in h2_d if d % 4 == 3])
    count2 = len(d_type2_h1) + len(d_type2_h2)
    print(f"Found {count2} HFDs of this type (all are maximal orders).")
    print(f"  d values (h=1): {d_type2_h1}")
    print(f"  d values (h=2): {d_type2_h2}")
    print(f"\nTotal HFDs of Type 2 = {len(d_type2_h1)} + {len(d_type2_h2)} = {count2}\n")

    # --- Final Calculation ---
    total_count = count1 + count2
    print("--- Total Count ---")
    print("The total number of such rings is the sum of the counts for the two disjoint types.")
    print(f"Final Equation: {count1} + {count2} = {total_count}")
    
    return total_count

# Execute the function to get the answer.
final_answer = solve_unique_length_factorizations()
# The final answer is printed by the function, but we also return it.
# print(f"\nFinal Answer: {final_answer}")
# The required format is <<<answer>>>
print(f"\n<<<{final_answer}>>>")
