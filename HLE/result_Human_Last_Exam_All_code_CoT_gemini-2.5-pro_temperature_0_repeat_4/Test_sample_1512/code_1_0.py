import math

# Wuxing architecture specifications
INT_SIZE_D = 5  # An 'int' is 5 decimal digits (D)

# Problem specifications
P_MAX_DIGITS = 100
Q_MAX_DIGITS = 100
# The product of two N-digit numbers can have at most 2N digits.
O_MAX_DIGITS = P_MAX_DIGITS + Q_MAX_DIGITS

# --- Calculations for minimized memory usage ---

# 1. Memory for p
# To store a 100-digit number, we can use an array of 5-digit integers.
p_chunks = math.ceil(P_MAX_DIGITS / INT_SIZE_D)
mem_p = p_chunks * INT_SIZE_D

# 2. Memory for q
# Same calculation for q.
q_chunks = math.ceil(Q_MAX_DIGITS / INT_SIZE_D)
mem_q = q_chunks * INT_SIZE_D

# 3. Memory for o (the product)
# The result can have up to 200 digits.
o_chunks = math.ceil(O_MAX_DIGITS / INT_SIZE_D)
mem_o = o_chunks * INT_SIZE_D

# 4. Total minimized memory m
total_mem_m = mem_p + mem_q + mem_o

# --- Output the results as an equation ---
# The prompt requires printing each number in the final equation.
print(f"Minimized memory for p ({P_MAX_DIGITS} digits) = {mem_p}D")
print(f"Minimized memory for q ({Q_MAX_DIGITS} digits) = {mem_q}D")
print(f"Minimized memory for o ({O_MAX_DIGITS} digits) = {mem_o}D")
print(f"Total minimized memory m = {mem_p} + {mem_q} + {mem_o} = {total_mem_m}D")