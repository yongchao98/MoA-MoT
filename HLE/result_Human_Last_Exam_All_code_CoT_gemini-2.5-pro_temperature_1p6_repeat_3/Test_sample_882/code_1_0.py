import math

def solve_and_calculate():
    """
    This function solves the game theory problem to find the value of p
    and then calculates the final requested expression.
    """
    # Parameters from the problem
    N = 8

    # The equation to solve for p is: p = 1 - (1 - p/N)^(3*N)
    # We reformulate it as finding the root of f(p) = 0
    # f(p) = 1 - (1 - p/N)^(3*N) - p
    def f(p):
        return 1.0 - (1.0 - p / N)**(3.0 * N) - p

    # We use the bisection method to find the non-trivial root of f(p) = 0 in (0, 1).
    # We choose a starting interval [low, high] where we know the root lies.
    # f(0.1) is positive and f(1.0) is negative, so a root is guaranteed in between.
    low = 0.1
    high = 1.0
    
    # Iterate 100 times for high precision, sufficient for this problem
    for _ in range(100):
        mid = (low + high) / 2.0
        if f(mid) > 0:
            low = mid
        else:
            high = mid
            
    # The value of p is the midpoint of the final interval
    p = (low + high) / 2.0
    
    # Calculate the expression requested by the user
    one_minus_p = 1 - p
    value = 10000 * one_minus_p
    final_answer = math.floor(value)
    
    # Output the required numbers from the final equation step-by-step
    print(f"For N = {N}, we solve the equation: p = 1 - (1 - p/{N})**(3*{N})")
    print(f"The calculated probability p is: {p:.6f}")
    print(f"The value of (1-p) is: {one_minus_p:.6f}")
    print(f"The value of 10000 * (1-p) is: {value:.4f}")
    print(f"The final result is floor({value:.4f}) = {final_answer}")

# Execute the function to print the solution
solve_and_calculate()