import sys

def calculate_c_k_formula(n, k):
    """
    Calculates the exact Shapley value c_k for player p_k using the derived formula.

    The plan is as follows:
    1. Define helper functions for the sum of the first n integers (S1) and the sum of the first n squares (S2).
    2. Implement the derived formula for the Shapley value c_k:
       c_k = k * S1(n)^3 - k^2 * S1(n)^2 + k * S1(n) * S2(n)
    3. Print the intermediate values and the final result as requested.
    """
    # Helper functions for sums of powers
    def sum_powers_1(num):
        """Calculates S1(n) = 1 + 2 + ... + n"""
        return num * (num + 1) // 2

    def sum_powers_2(num):
        """Calculates S2(n) = 1^2 + 2^2 + ... + n^2"""
        return num * (num + 1) * (2 * num + 1) // 6

    # Step 1: Calculate S1(n) and S2(n)
    s1_n = sum_powers_1(n)
    s2_n = sum_powers_2(n)

    # Step 2: Calculate the terms of the formula for c_k
    term1 = k * (s1_n ** 3)
    term2 = (k ** 2) * (s1_n ** 2)
    term3 = k * s1_n * s2_n

    c_k = term1 - term2 + term3

    # Step 3: Print the details of the calculation
    print(f"Calculating the fair division c_k for player p_{k}:")
    print("-" * 40)
    print(f"Total number of people (n): {n}")
    print(f"Player index (k): {k}")
    print("\nIntermediate values based on n:")
    print(f"  S1(n) = 1 + ... + n = {s1_n}")
    print(f"  S2(n) = 1^2 + ... + n^2 = {s2_n}")
    print("\nApplying the formula for c_k:")
    print(f"  c_k = k * S1(n)^3 - k^2 * S1(n)^2 + k * S1(n) * S2(n)")
    # "output each number in the final equation"
    print(f"  c_{k} = {k} * {s1_n}^3 - {k}^2 * {s1_n}^2 + {k} * {s1_n} * {s2_n}")
    print(f"  c_{k} = {term1} - {term2} + {term3}")
    print("\nFinal Result:")
    print(f"  The exact amount for player p_{k} is: ${c_k}")
    print("-" * 40)


if __name__ == '__main__':
    # Use command-line arguments if provided, otherwise use example values.
    if len(sys.argv) == 3:
        try:
            input_n = int(sys.argv[1])
            input_k = int(sys.argv[2])
            if not (input_n > 1 and 1 <= input_k <= input_n):
                print("Error: Ensure n > 1 and 1 <= k <= n.", file=sys.stderr)
                sys.exit(1)
            n_val, k_val = input_n, input_k
        except ValueError:
            print("Error: n and k must be integers.", file=sys.stderr)
            sys.exit(1)
    else:
        # Example values
        print("Usage: python your_script_name.py <n> <k>")
        print("Using example values n=4, k=2\n")
        n_val, k_val = 4, 2
        
    calculate_c_k_formula(n_val, k_val)