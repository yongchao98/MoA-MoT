import math

def solve_for_p():
    """
    Solves the equation 1 - (1 - p/N)^(3N) = p for p, given N=8.
    """
    N = 8
    
    # We need to solve f(p) = 1 - (1 - p/N)^(3*N) - p = 0
    # This is equivalent to p + (1 - p/N)^(3*N) - 1 = 0
    def f(p, n_val):
        return 1 - (1 - p / n_val)**(3 * n_val) - p

    # Use binary search to find the root for p in the interval (0, 1)
    low = 0.0
    high = 1.0
    # We iterate for a fixed number of steps to ensure high precision
    for _ in range(100):
        mid = (low + high) / 2
        if f(mid, N) > 0:
            low = mid
        else:
            high = mid
    
    p_solution = (low + high) / 2
    return p_solution

# 1. Solve for p
p = solve_for_p()
N_val = 8
three_N = 3 * N_val

# 2. Print the equation and verification
print("The problem reduces to solving the equation for p:")
print(f"1 - (1 - p / {N_val})**{three_N} = p")
print("\nNumerically solving this equation yields:")
# We need to print each number in the final equation. 
# Here we verify the equation with our solution for p.
lhs = 1 - (1 - p / N_val)**three_N
print(f"p = {p:.10f}")
print(f"LHS = 1 - (1 - {p:.10f} / {N_val})**{three_N} = {lhs:.10f}")
print(f"RHS = {p:.10f}")

# 3. Calculate the final value
val_1_minus_p = 1 - p
multiplied_val = 10000 * val_1_minus_p
final_answer = math.floor(multiplied_val)

print("\nCalculating the final answer:")
print(f"The value of (1-p) is: {val_1_minus_p:.10f}")
print(f"The value of 10000 * (1-p) is: {multiplied_val:.10f}")
print(f"The floor of this value is: {final_answer}")
