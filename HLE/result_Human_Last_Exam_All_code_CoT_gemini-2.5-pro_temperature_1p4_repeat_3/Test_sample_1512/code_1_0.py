import sys

# Step 1: Analyze memory requirements for p, q, and o on the Wuxing architecture.
# The Wuxing architecture uses a decimal system. The available data types relevant for storage are:
# - char: 1 decimal digit (1D)
# - int:  5 decimal digits (5D)
#
# The problem is to multiply two integers, p and q, each with a maximum of 100 digits.

# Step 2: Calculate the minimum memory for each number.
# The memory required to store a number is determined by its maximum number of digits.
# We must allocate memory for the worst-case scenario.

# Memory for p:
# p has a maximum of 100 digits. Storing 100 decimal digits requires 100 D of memory.
# This could be implemented as a 'char p[100];' in the specialized C.
mem_p = 100

# Memory for q:
# q also has a maximum of 100 digits, requiring 100 D of memory.
# This could be implemented as a 'char q[100];'.
mem_q = 100

# Memory for o (the output):
# The product of two 100-digit numbers can have up to 100 + 100 = 200 digits.
# Storing the result 'o' requires 200 D of memory.
# This could be implemented as a 'char o[200];'.
mem_o = 200

# Step 3: Calculate the total minimized memory 'm'.
# 'm' is the total memory use for p, q, and o. This is the sum of their
# individual memory requirements, assuming all three must reside in memory.
m = mem_p + mem_q + mem_o

# Step 4: Output the final equation as requested.
# The prompt asks to "output each number in the final equation".
# Here is the equation for the total memory calculation.
print(f"{mem_p} + {mem_q} + {mem_o} = {m}")

# Note: The problem mentions an "optimized C program". An optimized implementation
# would likely use arrays of 'int' (e.g., 'int p_chunks[20]') to perform calculations
# faster. However, the total storage in decimal digits remains the same:
# 20 chunks * 5 D/int = 100 D.
# Our calculation of m = 400 D is based on the fundamental data storage needs
# and is thus the minimized value.
