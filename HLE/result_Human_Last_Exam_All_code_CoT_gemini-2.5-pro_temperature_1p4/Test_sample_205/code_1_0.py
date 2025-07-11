import math

# Step 1: Analyze the memory requirements for calculating 100! on the WXVM.
# The result of 100! has 158 decimal digits. We must store it in an array.
# Using 'unsigned int' (4D) for each array element is efficient, as each can hold a 4-digit number (0-9999).
# - Size of result array: ceil(158 digits / 4 digits_per_int) = 40 elements.
# - Memory for result array: 40 * 4D = 160D.
#
# We also need variables for the calculation loop:
# - Outer loop counter 'i' (up to 100): a 'short' (2D) is sufficient.
# - Inner loop counter 'j' (up to 39): a 'short' (2D) is sufficient.
# - Array size tracker 'size' (up to 40): a 'short' (2D) is sufficient.
# - Carry variable for multiplication (up to ~1,000,000): needs a 'long' (8D).
#
# Step 2: Calculate the total minimum memory for variables.
# Total memory (z) = 160 (result) + 2 (i) + 2 (j) + 2 (size) + 8 (carry) = 174D.
z = 174

# Step 3: Calculate the value of 100!
# Python's math library can handle this large number calculation directly.
C = math.factorial(100)

# Step 4: Format the final output as z[C].
# The problem asks for the answer in this specific format.
print(f"{z}[{C}]")
