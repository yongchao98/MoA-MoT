from sympy import symbols, Eq

# Define the variables
A, B, C, D, E, F, G, H, I, J, K = symbols('A B C D E F G H I J K')

# Define the equations
equations = [
    Eq(F - K, 36),
    Eq(B - K, 111),
    Eq(B + E, 230),
    Eq(H + J, 25),
    Eq(E - A, 77),
    Eq(H + I, 36),
    Eq(C + E, 305),
    Eq(C - I, 205),
    Eq(C + J, 234),
    Eq(G, 3 * H),
    Eq(A + K, 42)
]

# The set of numbers
numbers = [3, 9, 16, 20, 39, 48, 75, 80, 121, 150, 225]

# Function to check if current assignment satisfies all equations
def satisfies_all(equations, assignment):
    return all(eq.subs(assignment) for eq in equations)

# Function to order variables by the number of constraints they are involved in
def order_variables_by_constraints(variables, equations):
    constraint_count = {var: 0 for var in variables}
    for eq in equations:
        for var in eq.free_symbols:
            if var in constraint_count:
                constraint_count[var] += 1
    return sorted(variables, key=lambda var: -constraint_count[var])

# Backtracking function with forward checking
def backtrack(assignment, variables, numbers, equations):
    if len(assignment) == len(variables):
        if satisfies_all(equations, assignment):
            return assignment
        return None

    # Get the next variable to assign
    var = variables[len(assignment)]

    for num in numbers:
        if num not in assignment.values():
            # Try assigning the number to the variable
            assignment[var] = num

            # Forward checking: check if this assignment is consistent
            if satisfies_all(equations, assignment):
                result = backtrack(assignment, variables, numbers, equations)
                if result is not None:
                    return result

            # Backtrack
            del assignment[var]

    return None

# Order variables by constraints
ordered_variables = order_variables_by_constraints([A, B, C, D, E, F, G, H, I, J, K], equations)

# Solve the problem using backtracking with forward checking
assignment = {}
solution = backtrack(assignment, ordered_variables, numbers, equations)

# If a solution is found, print it
if solution:
    result = [solution[A], solution[B], solution[C], solution[D],
              solution[E], solution[F], solution[G], solution[H],
              solution[I], solution[J], solution[K]]
    print(f"<<<{result}>>>")
else:
    print("No valid solution found.")