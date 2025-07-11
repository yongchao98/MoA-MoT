import math

def solve_quadratic_and_find_d(case_name, a, b, c):
    """Solves a quadratic equation and finds the root that satisfies d > 1."""
    print(f"Analyzing {case_name}...")
    # The prompt requires printing the numbers in the equation
    # To handle negative coefficients nicely, we format the string
    b_sign = "-" if b < 0 else "+"
    c_sign = "-" if c < 0 else "+"
    print(f"The derived equation is: {a} * d^2 {b_sign} {abs(b)} * d {c_sign} {abs(c)} = 0\n")

    # Calculate the discriminant
    discriminant = b**2 - 4 * a * c
    
    if discriminant < 0:
        print("No real solutions for d exist.\n")
        return None
    else:
        # Calculate the two possible roots for d
        d1 = (-b + math.sqrt(discriminant)) / (2 * a)
        d2 = (-b - math.sqrt(discriminant)) / (2 * a)
        print(f"The possible solutions for d are {d1:.2f} and {d2:.2f}")
        
        # Check which solution satisfies d > 1
        valid_d = None
        if d1 > 1:
            print(f"Solution d = {d1:.2f} is valid because it is greater than 1.")
            valid_d = d1
        elif d2 > 1:
            print(f"Solution d = {d2:.2f} is valid because it is greater than 1.")
            valid_d = d2
        else:
            print("Neither solution is greater than 1.")
        
        print("-" * 30)
        return valid_d

# --- Main execution ---
final_d = None

# Case 1: a1 = d, leading to 50d^2 - d - 51 = 0
d_case1 = solve_quadratic_and_find_d("Case 1 (a1 = d)", 50, -1, -51)
if d_case1 is not None:
    final_d = d_case1

# Case 2: a1 = 2d, leading to 51d^2 - d - 50 = 0
d_case2 = solve_quadratic_and_find_d("Case 2 (a1 = 2d)", 51, -1, -50)
# Although we already have a valid solution, we complete the analysis.
if d_case2 is not None and final_d is None:
    final_d = d_case2

# Print the final conclusion
if final_d is not None:
    print(f"\nThe only valid case gives d = {final_d:.2f}")
else:
    print("\nNo solution was found that satisfies d > 1.")
