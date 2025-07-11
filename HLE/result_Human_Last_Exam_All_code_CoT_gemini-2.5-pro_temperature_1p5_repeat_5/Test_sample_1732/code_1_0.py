import cmath

def solve_and_check(a, b, c):
    """
    Solves the quadratic equation ad^2 + bd + c = 0 and filters for d > 1.
    """
    print(f"Solving the equation: {a}d^2 + {b}d + {c} = 0")
    
    # Calculate the discriminant
    delta = (b**2) - 4*(a*c)
    
    # Find two solutions
    sol1 = (-b - cmath.sqrt(delta))/(2*a)
    sol2 = (-b + cmath.sqrt(delta))/(2*a)
    
    valid_solutions = []
    
    # Check if solutions are real and greater than 1
    if sol1.imag == 0:
        sol1_real = sol1.real
        if sol1_real > 1:
            valid_solutions.append(sol1_real)
            
    if sol2.imag == 0:
        sol2_real = sol2.real
        if sol2_real > 1:
            valid_solutions.append(sol2_real)
            
    return valid_solutions

# --- Main execution ---
final_d = None

print("--- Case 1: a_1 = d ---")
# Equation is 50d^2 - d - 51 = 0
case1_solutions = solve_and_check(a=50, b=-1, c=-51)
if case1_solutions:
    print(f"Valid solution(s) from Case 1: {case1_solutions}")
    # In this problem, we expect a unique solution.
    final_d = case1_solutions[0]
else:
    print("No solutions from Case 1 satisfy d > 1.")

print("\n--- Case 2: a_1 = 2d ---")
# Equation is 51d^2 - d - 50 = 0
case2_solutions = solve_and_check(a=51, b=-1, c=-50)
if case2_solutions:
    print(f"Valid solution(s) from Case 2: {case2_solutions}")
    # If Case 1 didn't yield a result, this one would be the answer.
    if final_d is None:
        final_d = case2_solutions[0]
else:
    print("No solutions from Case 2 satisfy d > 1.")

if final_d is not None:
    print(f"\nThe final value of d is: {final_d}")
else:
    print("\nNo valid solution for d was found.")
