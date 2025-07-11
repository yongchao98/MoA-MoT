def solve_trials():
    """
    Calculates and prints the minimum number of trials T(n) = 2n - 1 for given values of n.
    """
    
    n_values = [2, 3, 1234, 6712]
    results = []
    
    print("Based on the formula T(n) = 2n - 1:")
    
    for n in n_values:
        # According to the derived formula, T(n) = 2 * n - 1
        result = 2 * n - 1
        results.append(result)
        # The problem asks to output each number in the final equation
        print(f"T({n}) = 2 * {n} - 1 = {result}")
        
    # The final answer should be the comma-separated values
    final_answer = ", ".join(map(str, results))
    print(f"\nThe requested values T(2), T(3), T(1234), T(6712) are: {final_answer}")

solve_trials()