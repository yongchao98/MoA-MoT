from sympy import symbols, Eq, solve

# Define the symbols
A, C, D, E, F, G, H, I = symbols('A C D E F G H I')

# Define the equations based on the constraints
eq1 = Eq(D + E, 15)
eq2 = Eq(A + F, 201)
eq3 = Eq(D, 1.5 * E)
eq4 = Eq(C, 4.0 * D)
eq5 = Eq(E - I, -26)
eq6 = Eq(E, 2.0 * G)
eq7 = Eq(C + I, 68)

# Solve the equations
solution = solve((eq1, eq2, eq3, eq4, eq5, eq6, eq7), (A, C, D, E, F, G, I))

# Extract the values for the solved letters
solved_values = {key: solution[key] for key in solution}

# The numbers to assign
numbers = [3, 6, 9, 20, 32, 36, 80, 120, 121]

# Assign the solved values to the letters
assigned_values = {letter: solved_values[letter] for letter in solved_values}

# Remove the assigned numbers from the list
remaining_numbers = [num for num in numbers if num not in assigned_values.values()]

# Assign the remaining numbers to the remaining letters (B and H)
assigned_values['B'] = remaining_numbers[0]
assigned_values['H'] = remaining_numbers[1]

# Sort the values in alphabetical order
values = [assigned_values[letter] for letter in sorted(assigned_values.keys())]
print(values)