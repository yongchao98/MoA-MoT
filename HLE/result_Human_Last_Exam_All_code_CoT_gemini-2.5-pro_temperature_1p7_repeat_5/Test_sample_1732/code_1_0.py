import math

def solve_quadratic(a, b, c):
    """
    Solves the quadratic equation ax^2 + bx + c = 0 and returns real solutions.
    """
    # Calculate the discriminant
    delta = b**2 - 4*a*c
    
    # Check if real solutions exist
    if delta < 0:
        return []
        
    # Find two solutions
    sol1 = (-b - math.sqrt(delta)) / (2*a)
    sol2 = (-b + math.sqrt(delta)) / (2*a)
    
    return [sol1, sol2]

print("Based on the derivation, we have two possible cases for the value of d.")
print("We need to solve the quadratic equation for each case and check if the solution satisfies d > 1.")
print("-" * 50)

# --- Case 1 ---
print("Case 1: a_1 = d")
print("This case leads to the quadratic equation: 50d^2 - d - 51 = 0")
a1, b1, c1 = 50, -1, -51
print(f"Solving the equation where a={a1}, b={b1}, c={c1}...")
solutions_case1 = solve_quadratic(a1, b1, c1)
print(f"The mathematical solutions for d are: {solutions_case1}")

final_d = None
# Check for a valid solution in Case 1
for d_val in solutions_case1:
    if d_val > 1:
        final_d = d_val
        print(f"  -> Found a valid solution: d = {d_val} (since it is greater than 1)")
    else:
        print(f"  -> The solution d = {d_val} is not greater than 1.")

print("-" * 50)

# --- Case 2 ---
print("Case 2: a_1 = 2d")
print("This case leads to the quadratic equation: 51d^2 - d - 50 = 0")
a2, b2, c2 = 51, -1, -50
print(f"Solving the equation where a={a2}, b={b2}, c={c2}...")
solutions_case2 = solve_quadratic(a2, b2, c2)
print(f"The mathematical solutions for d are: {solutions_case2}")

# Check for a valid solution in Case 2
for d_val in solutions_case2:
    if d_val > 1:
        # This branch should not be hit based on our manual calculation, but we check anyway.
        final_d = d_val
        print(f"  -> Found a valid solution: d = {d_val} (since it is greater than 1)")
    else:
        print(f"  -> The solution d = {d_val} is not greater than 1.")

print("-" * 50)

# Final conclusion
if final_d is not None:
    print(f"\nThe only solution that satisfies the condition d > 1 is {final_d}.")
else:
    print("\nNo solution was found that satisfies all the problem's conditions.")
