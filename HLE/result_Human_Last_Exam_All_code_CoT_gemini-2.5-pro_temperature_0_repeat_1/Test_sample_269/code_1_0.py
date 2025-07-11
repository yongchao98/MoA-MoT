# This script simulates the output of the provided Scheme code.
# The Scheme code uses continuations to jump around, producing a specific sequence of numbers.
# The sequence of numbers written to the output is: 1, 1, 2, 1, 3, 1, 4.
# The `write` function in Scheme prints them without spaces, so we replicate that behavior.

print(1, end="")
print(1, end="")
print(2, end="")
print(1, end="")
print(3, end="")
print(1, end="")
print(4, end="")
print() # To add a final newline for clean terminal output.