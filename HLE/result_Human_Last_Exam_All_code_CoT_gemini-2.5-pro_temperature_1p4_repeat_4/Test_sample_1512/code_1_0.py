# Plan:
# 1. Define the memory size in digits (D) for the first input integer, p.
#    p has a max of 100 digits, so it requires 100 * 1D = 100D.
p_mem = 100

# 2. Define the memory size for the second input integer, q.
#    q also has a max of 100 digits, requiring 100 * 1D = 100D.
q_mem = 100

# 3. Define the memory size for the output integer, o.
#    The product of two 100-digit numbers can have up to 100 + 100 = 200 digits.
#    This requires 200 * 1D = 200D.
o_mem = 200

# 4. Calculate the total minimized memory use, m.
m = p_mem + q_mem + o_mem

# 5. Print the final equation showing each component as requested.
# The final equation for the total memory m is p_mem + q_mem + o_mem = m.
print(f"{p_mem} + {q_mem} + {o_mem} = {m}")