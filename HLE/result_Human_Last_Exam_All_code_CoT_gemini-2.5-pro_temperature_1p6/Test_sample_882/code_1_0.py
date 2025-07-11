import math

def solve_for_p():
    """
    Solves for the probability p in the mixed strategy Nash Equilibrium
    and computes the final requested value.
    """
    N = 8
    TOTAL_PLAYERS = 3 * N
    
    # The equilibrium equation is 1 - p = (1 - p/N)^(3N).
    # Let x = 1 - p. The equation becomes x = (1 - (1-x)/N)^(3N)
    # x = ((N - 1 + x)/N)^(3N)
    # For N=8, x = ((7+x)/8)^24

    print(f"The number of races is N = {N}")
    print(f"The number of players is 3N = {TOTAL_PLAYERS}")
    print("The equilibrium condition is given by the equation: 1 - p = (1 - p/N)^(3*N)")
    print("Substituting N = 8, we get: 1 - p = (1 - p/8)^24")
    print("Let x = 1 - p. The equation to solve for x is: x = ((7 + x) / 8) ^ 24")
    
    # We solve for x using fixed-point iteration: x_new = g(x_old)
    # Start with an initial guess for x. Since p is close to 1, x will be small.
    # Let's start with x = 0.
    x = 0.0
    
    # Iterate to find the fixed point with high precision.
    # 20 iterations are more than enough for convergence.
    for i in range(20):
        x_new = ((7 + x) / 8) ** 24
        # A simple check for convergence, though not strictly necessary here.
        if abs(x_new - x) < 1e-15:
            x = x_new
            break
        x = x_new

    one_minus_p = x
    p = 1 - one_minus_p

    print(f"\nSolved value for p (to six significant digits): {p:.6f}")
    print(f"Solved value for (1-p): {one_minus_p}")

    # Calculate the final result
    result = math.floor(10000 * one_minus_p)
    
    print(f"\nThe calculation is floor(10000 * (1-p)) = floor(10000 * {one_minus_p})")
    print(f"Result = floor({10000 * one_minus_p}) = {result}")
    
    return result

# Run the solver and print the final answer
final_answer = solve_for_p()
# The final result in the requested format will be added after this block
# For now, let's just make sure the print is clean.
# print(f"\nFinal Answer: {final_answer}")
print(f"\n<<< {final_answer} >>>")
