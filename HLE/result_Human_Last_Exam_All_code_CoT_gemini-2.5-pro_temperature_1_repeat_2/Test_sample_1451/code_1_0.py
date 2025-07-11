def solve_weighing_puzzle():
    """
    Calculates T(n) for a given list of n values based on the derived formulas.
    T(1) = 1
    T(n) = 2n - 1 for even n >= 2
    T(n) = 2n - 2 for odd n >= 3
    """
    
    n_values = [2, 3, 1234, 6712]
    results = []
    
    print("Calculating T(n) for each value based on the derived formula:")
    
    for n in n_values:
        if n == 1:
            # Handle the base case n=1
            result = 1
            print(f"T({n}) is a base case, and the result is {result}.")
            results.append(result)
        elif n % 2 == 0:
            # Formula for even n
            result = 2 * n - 1
            print(f"For even n={n}, the equation is T({n}) = 2 * {n} - 1 = {result}.")
            results.append(result)
        else: # n is odd and n > 1
            # Formula for odd n
            result = 2 * n - 2
            print(f"For odd n={n}, the equation is T({n}) = 2 * {n} - 2 = {result}.")
            results.append(result)
            
    # The final answer is the comma-separated list of results.
    final_answer_string = ",".join(map(str, results))
    print("\nThe final result is:")
    print(final_answer_string)

solve_weighing_puzzle()