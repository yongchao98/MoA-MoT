import math

def solve_cubic_equation():
    """
    Solves the equation p^3 + p - 1 = 0 for p in (0, 1) using the bisection method.
    """
    # The function whose root we want to find
    f = lambda p: p**3 + p - 1
    
    # Bisection method initialization
    low = 0.0
    high = 1.0
    
    # Iterate to find the root with high precision
    for _ in range(100):
        mid = (low + high) / 2
        if f(mid) < 0:
            low = mid
        else:
            high = mid
            
    # The root of the equation
    return (low + high) / 2

# 1. Solve for the optimal probability p0
p0 = solve_cubic_equation()

# 2. Calculate the constant C = 1 / (-ln(p0))
# This constant is the coefficient for the n*ln(n) term in the cost function.
constant = 1 / (-math.log(p0))

# 3. Print the final answer rounded to 3 decimal places.
print(f"The equation to optimize the questioning strategy is p^3 + p - 1 = 0.")
print(f"The optimal probability for a 'yes' answer is p0 = {p0:.6f}.")
print(f"The minimal cost is asymptotically C * n * ln(n), where C = 1 / (-ln(p0)).")
print(f"The value of the constant C is approximately {constant:.6f}.")
print("\nThe final answer, rounded to 3 decimal places, is:")
print(f"{constant:.3f}")
