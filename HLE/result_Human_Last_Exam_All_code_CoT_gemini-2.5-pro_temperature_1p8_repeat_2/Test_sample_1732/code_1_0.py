import math

def solve_quadratic_and_filter(a, b, c):
    """Solves the quadratic equation ax^2 + bx + c = 0 and filters for solutions > 1."""
    
    # Print the equation based on the prompt's requirement
    print(f"Solving the quadratic equation: ({a})d^2 + ({b})d + ({c}) = 0")
    
    # Calculate the discriminant
    delta = b**2 - 4*a*c
    
    if delta < 0:
        print("No real solutions.")
        return []
    
    # Calculate the two roots
    d1 = (-b + math.sqrt(delta)) / (2*a)
    d2 = (-b - math.sqrt(delta)) / (2*a)
    
    print(f"The mathematical solutions for d are: {d1} and {d2}")
    
    # Filter based on the condition d > 1
    valid_solutions = []
    if d1 > 1:
        valid_solutions.append(d1)
    if d2 > 1:
        valid_solutions.append(d2)
        
    if not valid_solutions:
        print("No solution found that satisfies the condition d > 1.")
    else:
        print(f"The valid solution(s) (d > 1) is/are: {valid_solutions}")
        
    return valid_solutions

def find_d():
    """
    Solves the problem by analyzing two cases derived from the problem's conditions.
    """
    final_d = None

    print("--- Analyzing Case 1 (a1 = d) ---")
    # This case leads to the equation 50d^2 - d - 51 = 0.
    solutions_case1 = solve_quadratic_and_filter(50, -1, -51)
    if solutions_case1:
        final_d = solutions_case1[0]
        
    print("\n--- Analyzing Case 2 (a1 = 2d) ---")
    # This case leads to the equation 51d^2 - d - 50 = 0.
    solutions_case2 = solve_quadratic_and_filter(51, -1, -50)
    # This case does not produce a valid result for d > 1, so we don't need to assign it.
    
    print("\n--- Conclusion ---")
    if final_d is not None:
        print(f"The only valid value for d from the analysis is {final_d}.")
    else:
        print("No valid solution was found.")

# Run the solver
find_d()
