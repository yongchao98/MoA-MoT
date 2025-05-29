from sympy import symbols, And, Or, Not, Implies, Equivalent, satisfiable

# Define the symbols
matthew_is_saint, david_is_saint = symbols('matthew_is_saint david_is_saint')

# Matthew's statement: "Matthew is a saint if and only if David is a saint"
matthew_statement = Equivalent(matthew_is_saint, david_is_saint)

# David's statement: "If Matthew is a saint then Matthew is a sinner"
david_statement = Implies(matthew_is_saint, Not(matthew_is_saint))

# Check the satisfiability of the statements
solution = satisfiable(And(matthew_statement, david_statement), all_models=True)

# Convert the generator to a list and print the solutions
solutions_list = list(solution)
print(solutions_list)