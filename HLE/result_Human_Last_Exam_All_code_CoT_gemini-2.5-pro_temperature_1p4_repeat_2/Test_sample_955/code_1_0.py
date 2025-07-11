import math

def solve_for_n():
    """
    Solves the quadratic equation derived from the problem's conditions to find n.
    """
    # The derived quadratic equation is 4*n^2 - 175*n + 600 = 0
    a = 4
    b = -175
    c = 600

    # Calculate the discriminant
    discriminant = (b**2) - (4 * a * c)

    if discriminant < 0:
        print("The equation has no real solutions for n.")
        return

    # Calculate the two possible solutions for n using the quadratic formula
    sqrt_discriminant = math.sqrt(discriminant)
    sol1 = (-b - sqrt_discriminant) / (2 * a)
    sol2 = (-b + sqrt_discriminant) / (2 * a)

    print("The relationship between the grid size 'n' and the probability is captured by the quadratic equation:")
    # Here we output each number in the final equation as requested.
    print(f"{a}*n^2 + ({b})*n + {c} = 0")
    
    n_solution = None
    # Check which solution is a positive, even integer, as required by the problem.
    for sol in [sol1, sol2]:
        # Check if the solution is an integer and meets the conditions
        if sol > 0 and sol == int(sol) and int(sol) % 2 == 0:
            n_solution = int(sol)
            break

    if n_solution is not None:
        print(f"\nSolving the equation gives two possible values: {sol1:.2f} and {sol2:.2f}.")
        print(f"The value that fits the condition of being an even integer is {n_solution}.")
        print(f"\nTherefore, the value of n is {n_solution}.")
    else:
        print("\nNo solution matching the problem's criteria (positive, even integer) was found.")

solve_for_n()