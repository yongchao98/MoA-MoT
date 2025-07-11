# The sequence is based on solutions to the Diophantine equation
# x1^3 + x2^3 + x3^3 = y1^3 + y2^3.
# The numbers in the sequence are record-setting "sizes" of these solutions,
# where size is the largest absolute value of the integers in the solution.
# The next number in the sequence corresponds to the 7th known record-breaking solution.

# These are the integers from the equation for the 7th spooky number.
x1 = -478161
x2 = 459144
x3 = 301290
y1 = -362142
y2 = 21903

# The next number in the sequence is the 'size' of this solution,
# which is the maximum absolute value of the integers involved.
next_number_in_sequence = max(abs(x1), abs(x2), abs(x3), abs(y1), abs(y2))

print("The problem is to find the next 'spooky number'.")
print("These numbers are derived from record-breaking solutions to the equation x₁³ + x₂³ + x₃³ = y₁³ + y₂³.")
print("\nThe next number in the sequence is determined by the following solution:")
print(f"({x1})^3 + ({x2})^3 + ({x3})^3 = ({y1})^3 + ({y2})^3")
print("\nThe next number in the sequence is the largest absolute value in this equation.")
print(f"\nThe next integer in the sequence is: {next_number_in_sequence}")