import math

def solve_for_p():
    """
    Numerically solves the equilibrium equation to find p and calculates the final result.
    The equation to solve is p - (1 - (1 - p/8)**24) = 0.
    For cleaner syntax, we solve f(p) = p + (1 - p/8)**24 - 1 = 0.
    """
    
    # Define the function f(p) whose root we want to find.
    def f(p):
        N = 8
        NUM_PLAYERS = 3 * N
        return p + (1 - p / N)**(NUM_PLAYERS) - 1

    # Bisection method to find the non-trivial root of f(p) = 0 in (0, 1)
    low = 1e-9  # Start slightly above 0 to avoid the trivial root
    high = 1.0
    
    # f(low) is positive, f(high) is negative, so a root exists in between.
    for _ in range(100):  # 100 iterations are more than enough for high precision
        mid = (low + high) / 2
        if f(mid) * f(low) > 0:
            low = mid
        else:
            high = mid
            
    p = (low + high) / 2
    
    # Parameters for printing the final equation
    N = 8
    num_players = 3 * N
    
    lhs = p
    rhs = 1 - (1 - p / N)**num_players

    print("The equilibrium condition is that the payoff for the discrete strategy equals 1/3.")
    print(f"This leads to the equation: p = 1 - (1 - p/{N})^({num_players})")
    print(f"Solving numerically for p, we get p ≈ {p:.6f}")
    print(f"Plugging p back into the equation:")
    print(f"{lhs:.6f} ≈ 1 - (1 - {p:.6f}/{N})^{num_players}")
    print(f"{lhs:.6f} ≈ {rhs:.6f}")

    # Calculate the final value as requested
    result = math.floor(10000 * (1 - p))
    print("\nThe final calculation is floor(10000 * (1 - p))")
    print(f"p = {p}")
    print(f"1 - p = {1-p}")
    print(f"10000 * (1 - p) = {10000 * (1 - p)}")
    print(f"floor(10000 * (1 - p)) = {result}")
    
    return result

# Execute the solver and print the final answer
final_answer = solve_for_p()
# The final result will be printed by the function above.
# print(f"<<<{final_answer}>>>") # As per user instructions, output the answer at the end
# Since the user asked me not to copy-paste, printing directly is better.
# However, the final answer needs the specific format.

print(f"\nFinal Answer:")
print(f'<<<{final_answer}>>>')