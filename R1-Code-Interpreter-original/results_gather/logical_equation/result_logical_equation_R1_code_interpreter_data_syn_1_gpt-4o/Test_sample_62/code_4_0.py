from sympy import symbols, Eq, solve

# Define the symbols
A, B, C, D, E, G, H, I = symbols('A B C D E G H I')

# Define the equations
eq1 = Eq(C, 3.0 * E)
eq2 = Eq(E + I, 18)
eq3 = Eq(B + H, 12)
eq4 = Eq(B + C, 52)
eq5 = Eq(D - E, -13)
eq6 = Eq(B - A, -9)

# Solve the system of equations
solution = solve((eq1, eq2, eq3, eq4, eq5, eq6), (A, B, C, D, E, G, H, I))

# List of numbers to assign
numbers = [1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50]

# Extract the values in alphabetical order
values = [solution[A], solution[B], solution[C], solution[D], solution[E], None, None, solution[H], solution[I], None, None, None]

# Find the numbers for F, G, J, K, and L
used_numbers = set(values) - {None}
remaining_numbers = list(set(numbers) - used_numbers)

# Assign the remaining numbers to F, G, J, K, and L
values[5] = remaining_numbers[0]  # Assign one of the remaining numbers to F
values[6] = remaining_numbers[1]  # Assign another remaining number to G
values[9] = remaining_numbers[2]  # Assign another remaining number to J
values[10] = remaining_numbers[3]  # Assign another remaining number to K
values[11] = remaining_numbers[4]  # Assign the last remaining number to L

print(values)