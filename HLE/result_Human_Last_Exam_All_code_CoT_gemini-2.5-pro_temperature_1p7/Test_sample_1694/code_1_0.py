import sympy

# Define a function to pretty print the sum equation
def print_sum_equation(terms, result):
    """Formats and prints the sum equation."""
    equation_parts = []
    # Use pretty printing for each term
    first_term = sympy.pretty(terms[0])
    equation_parts.append(first_term)
    
    for term in terms[1:]:
        # Add parentheses for negative terms to make the sum clearer
        pretty_term = sympy.pretty(term)
        if pretty_term.startswith('-'):
            equation_parts.append(f"+ ({pretty_term})")
        else:
            equation_parts.append(f"+ {pretty_term}")
            
    print(" ".join(equation_parts), f"= {sympy.pretty(result)}")

# --- Equation 1 ---
# Solve for Y_11 in 11*Y_11 = -53/12
y11 = -sympy.Rational(53, 132)
# The first coordinates of X1 are the square roots of Y_11
x1_coord_sol1 = sympy.sqrt(y11)
x1_coord_sol2 = -sympy.sqrt(y11)

# --- Equation 2 ---
# Solve for Z_11 in 10*Z_11 = -3/11
z11 = -sympy.Rational(3, 110)
# The first coordinates of X2 are the square roots of Z_11
x2_coord_sol1 = sympy.sqrt(z11)
x2_coord_sol2 = -sympy.sqrt(z11)

# The "solutions" are the set of all possible matrices. We sum the first coordinate of each.
all_first_coords = [x1_coord_sol1, x1_coord_sol2, x2_coord_sol1, x2_coord_sol2]
total_sum = sum(all_first_coords)

# Print the final equation with each number explicitly shown.
print("The first coordinates of the solutions are:")
print(f"For X1: {sympy.pretty(x1_coord_sol1)} and {sympy.pretty(x1_coord_sol2)}")
print(f"For X2: {sympy.pretty(x2_coord_sol1)} and {sympy.pretty(x2_coord_sol2)}")
print("\nThe sum of the first coordinate of all solutions is:")
print_sum_equation(all_first_coords, total_sum)
