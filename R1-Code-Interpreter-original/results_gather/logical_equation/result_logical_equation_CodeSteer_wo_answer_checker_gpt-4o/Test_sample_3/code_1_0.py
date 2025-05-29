from sympy import symbols, Eq, nonlinsolve

# Define the variables
A, B, C, D, E, F, G, H, I = symbols('A B C D E F G H I', integer=True)

# Define the possible values
possible_values = {1, 3, 4, 9, 16, 27, 36, 80, 121}

# Define the equations
equations = [
    Eq(A - D, -77),
    Eq(I, 3 * E),
    Eq(E - H, -27),
    Eq(G - A, -2),
    Eq(H, 4 * E)
]

# Use nonlinsolve to find solutions
solutions = nonlinsolve(equations, (A, D, E, G, H, I))

# Filter solutions to satisfy inequalities and possible values
valid_solution = None
for sol in solutions:
    a, d, e, g, h, i = sol
    # Substitute the solution into inequalities
    if (h.subs({H: h, A: a}) > a) and (i.subs({I: i, B: B}) > B):
        # Check if the solution values are in the possible values set
        if all(value in possible_values for value in sol):
            # Assign remaining values
            remaining_values = possible_values - {a, d, e, g, h, i}
            # Assign B, C, F from remaining values
            for b in remaining_values:
                for c in remaining_values - {b}:
                    for f in remaining_values - {b, c}:
                        if i > b:
                            valid_solution = (a, b, c, d, e, f, g, h, i)
                            break
                    if valid_solution:
                        break
                if valid_solution:
                    break
    if valid_solution:
        break

# Output the solution in the required format
if valid_solution:
    result = list(valid_solution)
    print(f"<<<{result}>>>")
else:
    print("No valid solution found.")