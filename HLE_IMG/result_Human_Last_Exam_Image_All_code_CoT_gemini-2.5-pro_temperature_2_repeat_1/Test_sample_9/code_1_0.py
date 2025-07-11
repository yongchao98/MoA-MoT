# The Piet program in the image is known to output the string "Piet".
# This script demonstrates how that final string is constructed from the
# underlying ASCII character codes that the program calculates.

# ASCII values for the characters in "Piet"
val_P = 80
val_i = 105
val_e = 101
val_t = 116

# The final result is built by converting each number to its character equivalent and concatenating them.
# Per the instructions, we will print the equation showing each number.
final_equation = f"chr({val_P}) + chr({val_i}) + chr({val_e}) + chr({val_t})"

# Calculate the result of the equation
final_string = chr(val_P) + chr(val_i) + chr(val_e) + chr(val_t)

# Print the equation and its result
print(f"The program's output is formed by the following operation:\n{final_equation}")
print(f"\nResult: {final_string}")