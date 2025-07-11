# Plan: Calculate the minimum memory (m) in decimal digits (D) required to
# store two 100-digit integers (p, q) and their product (o).

# According to the Wuxing architecture, the 'char' type is 1D (1 decimal digit),
# which is the most memory-efficient way to store digits.

# 1. Determine the memory required for the first input integer, p.
# p has a maximum of 100 digits.
mem_p = 100  # in D

# 2. Determine the memory required for the second input integer, q.
# q has a maximum of 100 digits.
mem_q = 100  # in D

# 3. Determine the memory required for the output product, o.
# The product of two numbers with N and M digits has at most N + M digits.
# So, o has a maximum of 100 + 100 = 200 digits.
mem_o = mem_p + mem_q

# 4. Calculate the total minimized memory, m.
# This is the sum of the memory for p, q, and o.
total_memory = mem_p + mem_q + mem_o

# 5. Print the equation as requested.
# The 'f-string' formatting inserts the variable values into the string.
print(f"{mem_p} + {mem_q} + {mem_o} = {total_memory}")