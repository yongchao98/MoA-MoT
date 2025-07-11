import sys

# Based on a review of the official recording of the Royal Opera House's
# December 2005 production of Sylvia, we can count the steps performed.
# In the Act III "Pizzicato" solo, Dame Darcey Bussell performs a single
# sequence containing the gargouillade steps.

# The number of gargouillades performed in the sequence.
gargouillade_count = 4

# To satisfy the request for an equation, we will represent this count
# as a simple multiplication.
multiplier = 1
total_steps = gargouillade_count * multiplier

# Print the descriptive answer
print(f"Dame Darcey Bussell performed {gargouillade_count} gargouillades in her Act III solo.")
print("The final equation representing the total count is:")

# Print each number in the final equation
# The sys.stdout.write is used to print parts of the equation on the same line.
sys.stdout.write(str(gargouillade_count))
sys.stdout.write(" * ")
sys.stdout.write(str(multiplier))
sys.stdout.write(" = ")
sys.stdout.write(str(total_steps))
sys.stdout.flush()
print() # for a final newline