# The maximum number of digits for the input integer p.
# Using a char array is the most memory-efficient, so storage is 1D per digit.
p_storage_d = 100

# The maximum number of digits for the input integer q.
q_storage_d = 100

# The maximum number of digits for the product o = p * q is the sum of the
# number of digits in p and q.
o_storage_d = p_storage_d + q_storage_d

# The minimized total memory use 'm' is the sum of the storage needed for p, q, and o.
m_total_d = p_storage_d + q_storage_d + o_storage_d

# Print the final equation showing the storage for each number and the total.
print(f"{p_storage_d} + {q_storage_d} + {o_storage_d} = {m_total_d}")