import math

def solve_and_check(a, b, c, case_name):
    """Solves the quadratic equation ax^2 + bx + c = 0 and checks for valid solutions."""
    print(f"----- {case_name} -----")
    print(f"The equation for d is: {a}d^2 + ({b})d + ({c}) = 0")
    
    # Calculate the discriminant
    delta = b**2 - 4*a*c
    
    if delta < 0:
        print("The equation has no real roots.")
        return None
    
    # Calculate the two roots
    root1 = (-b + math.sqrt(delta)) / (2*a)
    root2 = (-b - math.sqrt(delta)) / (2*a)
    
    print(f"The roots are d = {root1} and d = {root2}")
    
    # Check the condition d > 1
    valid_solution = None
    if root1 > 1:
        print(f"Solution d = {root1} is valid (d > 1).")
        valid_solution = root1
    else:
        print(f"Solution d = {root1} is not valid (d is not > 1).")
        
    if root2 > 1:
        print(f"Solution d = {root2} is valid (d > 1).")
        valid_solution = root2
    else:
        print(f"Solution d = {root2} is not valid (d is not > 1).")
        
    return valid_solution

# --- Main ---
# Case 1: 50d^2 - d - 51 = 0
final_d = solve_and_check(50, -1, -51, "Case 1: a_1 = d")
print("\n")

# Case 2: 51d^2 - d - 50 = 0
solve_and_check(51, -1, -50, "Case 2: a_1 = 2d")
print("\n")

if final_d is not None:
    print("-----------------------------------")
    print(f"The only valid value for d that satisfies d > 1 is found in Case 1.")
    print(f"Final Answer: d = {final_d}")