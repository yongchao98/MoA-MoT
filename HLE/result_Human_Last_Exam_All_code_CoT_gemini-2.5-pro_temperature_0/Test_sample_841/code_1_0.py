import sys

def solve_and_print():
    """
    This function solves the problem by counting the number of specified rings
    that are Half-Factorial Domains (HFDs).
    """
    # Step 1: Define the known sets of d for which the class number h(-d) is 1 or 2.
    # These are fundamental results from algebraic number theory.
    h1_d = {1, 2, 3, 7, 11, 19, 43, 67, 163}
    h2_d = {5, 6, 10, 13, 15, 22, 35, 37, 51, 58, 91, 115, 123, 187, 235, 267, 403, 427}

    # --- Part 1: Count HFDs that are rings of integers O(Q(sqrt(-d))) ---
    # A ring of integers O(Q(sqrt(-d))) is an HFD if and only if h(-d) is 1 or 2.
    # The total number of such rings is the sum of the sizes of the two sets above.
    num_hfd_maximal_orders = len(h1_d) + len(h2_d)

    # --- Part 2: Count HFDs that are non-integrally closed rings Z[sqrt(-d)] ---
    # A ring Z[sqrt(-d)] is not integrally closed if d = 3 (mod 4).
    # Such a ring is an HFD if d=3 or h(-d)=2.

    # Case d=3:
    # d=3 satisfies d=3 (mod 4). This gives one HFD ring: Z[sqrt(-3)].
    count_d_is_3 = 1
    d_val_3 = 3

    # Case h(-d)=2:
    # We need to find which d from the h2_d set also satisfy d = 3 (mod 4).
    hfd_non_maximal_from_h2 = {d for d in h2_d if d % 4 == 3}
    count_h2_non_maximal = len(hfd_non_maximal_from_h2)

    # Total number of HFDs of this non-maximal type.
    num_hfd_non_maximal_orders = count_d_is_3 + count_h2_non_maximal

    # --- Final Calculation ---
    # The two sets of rings are disjoint:
    # Part 1 counts maximal orders O(Q(sqrt(-d))).
    # Part 2 counts non-maximal orders Z[sqrt(-d)] (since d=3 mod 4).
    # The total size is the sum of the two counts.
    total_size = num_hfd_maximal_orders + num_hfd_non_maximal_orders

    # --- Print the detailed explanation and result ---
    print("This program calculates the size of a set of number rings with unique factorization lengths.")
    print("\n--- Part 1: Counting HFDs that are Rings of Integers O(Q(sqrt(-d))) ---")
    print("A ring of integers is a Half-Factorial Domain (HFD) if its class number is 1 or 2.")
    print(f"The number of square-free d > 0 with class number 1 is {len(h1_d)}.")
    print(f"The number of square-free d > 0 with class number 2 is {len(h2_d)}.")
    print(f"Total number of HFDs of this type = {len(h1_d)} + {len(h2_d)} = {num_hfd_maximal_orders}")

    print("\n--- Part 2: Counting HFDs that are non-integrally closed rings Z[sqrt(-d)] ---")
    print("The ring Z[sqrt(-d)] is non-integrally closed when d = 3 (mod 4).")
    print("Such a ring is an HFD if d=3 or the class number h(-d)=2.")
    print(f"The case d={d_val_3} gives 1 HFD ring.")
    print(f"The number of d values with h(-d)=2 and d=3(mod 4) is {count_h2_non_maximal}.")
    print(f"Total number of HFDs of this type = {count_d_is_3} + {count_h2_non_maximal} = {num_hfd_non_maximal_orders}")

    print("\n--- Final Result ---")
    print("The two sets of HFDs are disjoint, so the total size is the sum of their sizes.")
    print(f"Total size = (Number of maximal HFDs) + (Number of non-maximal HFDs)")
    print(f"Total size = {num_hfd_maximal_orders} + {num_hfd_non_maximal_orders} = {total_size}")

    # The final answer in the requested format is printed outside the function
    # to ensure it's the last thing in the output stream.
    return total_size

if __name__ == '__main__':
    final_answer = solve_and_print()
    # The problem asks for the final answer in a specific format at the end.
    # To avoid printing it from a library, we'll use a simple print.
    sys.stdout.write(f"\n<<<{final_answer}>>>\n")
