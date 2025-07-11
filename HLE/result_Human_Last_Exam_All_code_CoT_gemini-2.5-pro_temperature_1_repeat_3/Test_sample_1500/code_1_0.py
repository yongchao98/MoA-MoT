def analyze_betti_numbers_for_sun(n):
    """
    This function demonstrates the counterexample for question (b)
    by comparing the second Betti number of a regular vs. a singular
    coadjoint orbit for the Lie group SU(n).
    """
    print(f"--- Analysis for G = SU({n}) ---")

    if n <= 2:
        print(f"For n={n}, the rank is {n-1}, and b2(CP^{n-1}) is 1. The values match.")
        print("Let's consider the more illustrative case where n > 2.")
        return

    # For a regular element lambda, the orbit is the full flag manifold,
    # and its second Betti number is the rank of the group.
    rank_sun = n - 1
    print(f"For a regular orbit, the expected second Betti number b2 is the rank of SU({n}), which is {n} - 1 = {rank_sun}.")

    # For a specific singular element lambda, the orbit is the complex projective space CP^{n-1}.
    # The second Betti number of CP^{n-1} is 1.
    b2_cp_n_minus_1 = 1
    print(f"For a singular orbit like the complex projective space CP^{n-1}, the second Betti number b2 is {b2_cp_n_minus_1}.")

    # The final equation check
    print("\nIs the second Betti number always n-1?")
    print(f"Comparing the two values: {rank_sun} and {b2_cp_n_minus_1}")
    if rank_sun != b2_cp_n_minus_1:
        print(f"Since {rank_sun} != {b2_cp_n_minus_1}, the statement is false for SU({n}).")
    else:
        # This case won't be hit due to the n > 2 check above.
        print(f"Since {rank_sun} == {b2_cp_n_minus_1}, the statement holds for SU({n}).")

def main():
    # Demonstrate the counterexample for n=3, which is a clear case.
    analyze_betti_numbers_for_sun(3)

    # Print the final answers in the requested format.
    print("\n--- Final Answers ---")
    final_answer = "(a) True; (b) No; (c) No"
    print(final_answer)
    print(f"\n<<<{final_answer}>>>")

if __name__ == "__main__":
    main()