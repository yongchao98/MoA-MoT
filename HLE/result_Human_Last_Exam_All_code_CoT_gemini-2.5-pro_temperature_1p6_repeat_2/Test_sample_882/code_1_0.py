import math

def solve_for_p():
    """
    Solves the equilibrium equation for p and computes the final result.
    """
    N = 8
    M = 3 * N

    # The equation to solve for p is: p = 1 - (1 - p/N)**M
    # We can solve this using fixed-point iteration p_k+1 = g(p_k)
    # where g(p) = 1 - (1 - p/N)**M
    
    # Start with an initial guess for p. The root is close to 1.
    p = 0.95
    
    # Iterate for a number of steps to ensure convergence to high precision.
    for _ in range(100):
        p_next = 1 - (1 - p / N)**M
        # Check if the value has converged
        if abs(p_next - p) < 1e-12:
            p = p_next
            break
        p = p_next

    # The problem asks to output each number in the final equation.
    # The final equation is the calculation of the result: floor(10000 * (1 - p)).
    # We will print the values used in this calculation.
    value_to_floor = 10000 * (1 - p)
    result = math.floor(value_to_floor)
    
    print(f"The probability p is: {p}")
    print(f"The value of 10000 * (1 - p) is: {value_to_floor}")
    print(f"The final integer result is: {result}")

solve_for_p()