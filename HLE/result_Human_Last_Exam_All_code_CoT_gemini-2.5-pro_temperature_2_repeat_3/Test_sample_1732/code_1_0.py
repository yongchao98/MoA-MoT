import math

def solve_quadratic(a, b, c):
    """Solves the quadratic equation ax^2 + bx + c = 0 and returns a list of real solutions."""
    delta = b**2 - 4*a*c
    if delta < 0:
        return []
    sqrt_delta = math.sqrt(delta)
    x1 = (-b + sqrt_delta) / (2*a)
    x2 = (-b - sqrt_delta) / (2*a)
    return [x1, x2]

# --- Case 1: a_1 = d ---
# Equation: 50d^2 - d - 51 = 0
case1_coeffs = (50, -1, -51)
solutions_case1 = solve_quadratic(*case1_coeffs)

# --- Case 2: a_1 = 2d ---
# Equation: 51d^2 - d - 50 = 0
case2_coeffs = (51, -1, -50)
solutions_case2 = solve_quadratic(*case2_coeffs)

# Check the solutions against the condition d > 1
final_d = None
final_coeffs = None

for d_val in solutions_case1:
    if d_val > 1:
        final_d = d_val
        final_coeffs = case1_coeffs

# If no valid solution was found in case 1, check case 2
if final_d is None:
    for d_val in solutions_case2:
        if d_val > 1:
            final_d = d_val
            final_coeffs = case2_coeffs

# Print the final results
if final_d is not None and final_coeffs is not None:
    print("From the problem conditions, we derived two possible quadratic equations for d.")
    print("The equation that yields a valid solution (where d > 1) is:")
    # Print out each number in the equation
    print(f"({final_coeffs[0]}) * d^2 + ({final_coeffs[1]}) * d + ({final_coeffs[2]}) = 0")
    print("\nSolving this equation gives two roots, but only one is greater than 1.")
    print(f"The value of d is: {final_d}")
else:
    print("No solution for d > 1 was found.")
