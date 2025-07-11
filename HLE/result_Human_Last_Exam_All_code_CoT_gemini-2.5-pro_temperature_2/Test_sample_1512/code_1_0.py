import math

# Define the maximum number of decimal digits (D) for the inputs p and q.
digits_p = 100
digits_q = 100

# The Wuxing C compiler provides a 5-digit 'int' type.
# This is the optimal chunk size for large number arithmetic on this architecture.
int_size_in_digits = 5

# The maximum number of digits in the output o = p * q is the sum of the digits of p and q.
digits_o = digits_p + digits_q

# Calculate the memory in D required to store each number.
# We represent each large number as an array of 5-digit 'int's.
# The memory required is the number of digits rounded up to the nearest multiple of our chunk size (5).
# Since 100 and 200 are perfect multiples of 5, no rounding is needed.
mem_p = math.ceil(digits_p / int_size_in_digits) * int_size_in_digits
mem_q = math.ceil(digits_q / int_size_in_digits) * int_size_in_digits
mem_o = math.ceil(digits_o / int_size_in_digits) * int_size_in_digits

# The minimized total memory use 'm' is the sum of the memory for p, q, and the result o.
m = mem_p + mem_q + mem_o

# Print the final calculation showing each component.
print(f"{int(mem_p)} + {int(mem_q)} + {int(mem_o)} = {int(m)}")