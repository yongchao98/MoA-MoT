# Define the maximum number of digits for the input integers p and q.
p_max_digits = 100
q_max_digits = 100

# In the Wuxing architecture, the 'char' data type is 1D (1 digit),
# which is the most memory-efficient way to store a single digit.
size_per_digit_in_D = 1

# Calculate the memory required for p.
mem_p = p_max_digits * size_per_digit_in_D

# Calculate the memory required for q.
mem_q = q_max_digits * size_per_digit_in_D

# The maximum number of digits in the product o = p * q is the sum
# of the number of digits in p and q.
o_max_digits = p_max_digits + q_max_digits

# Calculate the memory required for the product o.
mem_o = o_max_digits * size_per_digit_in_D

# Calculate the minimized total memory use 'm'.
m_total_mem = mem_p + mem_q + mem_o

# The final equation is m = mem_p + mem_q + mem_o.
# The instruction is to output each number in this final equation.
# We will print the equation in a clear format.
print(f"{mem_p} + {mem_q} + {mem_o} = {m_total_mem}")