# Wuxing architecture parameters
# The size of an 'int' data type is 5 decimal digits (5D).
int_size_in_digits = 5

# Problem parameters: maximum number of digits for inputs p and q.
p_max_digits = 100
q_max_digits = 100

# The optimal C program would represent the large numbers as arrays of 'int's.
# This aligns with the 5D register size for efficient processing.

# Calculate the memory required for p in decimal digits (D).
# An array of 100/5 = 20 integers is needed. Total memory = 20 * 5D.
p_memory_D = p_max_digits

# Calculate the memory required for q in decimal digits (D).
# An array of 100/5 = 20 integers is needed. Total memory = 20 * 5D.
q_memory_D = q_max_digits

# Calculate the memory required for the product o.
# The maximum number of digits in the product o is the sum of the digits of the operands.
o_max_digits = p_max_digits + q_max_digits
# An array of 200/5 = 40 integers is needed. Total memory = 40 * 5D.
o_memory_D = o_max_digits

# Calculate the minimized total memory use 'm' by summing the memory for p, q, and o.
m = p_memory_D + q_memory_D + o_memory_D

# Output the final equation showing how the total memory is calculated.
print(f"To store p, q, and their product o, we need:")
print(f"Memory for p: {p_memory_D} D")
print(f"Memory for q: {q_memory_D} D")
print(f"Memory for o: {o_memory_D} D")
print("\nThe minimized total memory use (m) is the sum:")
print(f"m = {p_memory_D} + {q_memory_D} + {o_memory_D} = {m}")