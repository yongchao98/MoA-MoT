import math

def solve_k_matching_complexity():
    """
    Determines the maximum integer k for which counting k-matchings can be done
    in subcubic time, based on fine-grained complexity assumptions.
    """
    # The current best known upper bound for the matrix multiplication exponent ω (omega).
    omega = 2.3728596

    print(
        "The time complexity for counting k-matchings is O(n^α(k)), where:\n"
        "α(k) = ω * floor(k/3) + (k mod 3)\n"
        "We need to find the maximum integer k such that α(k) < 3.\n"
        f"We will use ω ≈ {omega:.4f}.\n"
    )

    k = 1
    max_k_for_subcubic = 0

    while True:
        floor_k_div_3 = math.floor(k / 3)
        k_mod_3 = k % 3

        # Calculate the exponent α(k)
        alpha_k = omega * floor_k_div_3 + k_mod_3

        print(f"--- Checking k = {k} ---")
        print(f"The equation is: ω * floor({k}/3) + ({k} mod 3) < 3")
        print(f"Substituting values: {omega:.4f} * {floor_k_div_3} + {k_mod_3} < 3")
        print(f"Calculation: {alpha_k:.4f} < 3")

        if alpha_k < 3:
            print("Result: The condition is TRUE. Counting k-matchings is subcubic.")
            max_k_for_subcubic = k
            k += 1
        else:
            print("Result: The condition is FALSE. Counting k-matchings is not subcubic.")
            break

    print("\n" + "="*40)
    print(f"The maximum value of k for which the complexity is subcubic is {max_k_for_subcubic}.")
    print("="*40)


if __name__ == "__main__":
    solve_k_matching_complexity()
