import math

def explain_and_demonstrate_complexity():
    """
    Demonstrates the algorithm's complexity with a concrete example
    and prints the final analysis.
    """
    # 1. Define example inputs
    n = 137
    c = 10

    # 2. Trace the algorithm to count loop iterations
    # We don't need the full algorithm, just the part that determines iterations.
    iterations = 0
    if n > 0:
        q = n // c
        while q > 0:
            iterations += 1
            x = q
            q = x // c

    # 3. Perform the theoretical calculation for comparison
    # The number of while loop iterations is floor(log_c(n)) for n > 0
    log_val = math.log(n, c)
    theoretical_iterations = math.floor(log_val)

    # 4. Print the analysis and results
    print(f"Analyzing the complexity of IMAG(n, c) with example: n={n}, c={c}")
    print("-----------------------------------------------------------------")
    print("The algorithm's runtime is determined by the number of 'while' loop runs.")
    
    print(f"\nFor n={n} and c={c}:")
    print(f"  - The number of times the loop actually runs is: {iterations}")
    
    print(f"\nThe theoretical number of iterations is floor(log_c(n)). Let's calculate this:")
    # The following print statement shows the "final equation" with its numbers
    print(f"  - floor(log_{c}({n})) = floor({log_val:.4f}) = {theoretical_iterations}")

    print("\nResult:")
    print("The actual number of iterations matches the theoretical value. This demonstrates")
    print("that the number of operations is proportional to log_c(n).")
    
    print("\n-----------------------------------------------------------------")
    print("Final Conclusion:")
    print("The general computational time complexity of the algorithm is O(log_c(n)).")

# Run the demonstration
explain_and_demonstrate_complexity()