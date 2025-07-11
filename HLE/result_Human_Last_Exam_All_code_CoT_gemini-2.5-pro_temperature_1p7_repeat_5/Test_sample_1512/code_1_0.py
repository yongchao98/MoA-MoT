# The maximum number of digits for each input integer, p and q.
digits_p = 100
digits_q = 100

# The maximum number of digits for the output integer o, which is the product of p and q.
# The number of digits in a product is at most the sum of the number of digits of the operands.
digits_o = digits_p + digits_q

# In the Wuxing architecture, data types are measured in decimal digits (D).
# char = 1D, int = 5D, long = 10D.
# An optimized program needs to store the numbers. The total storage in D for a number
# is simply its number of digits, assuming a packed representation (like an array of chars
# or ints) with no wasted space.
# Therefore, the memory usage for each variable is equal to its maximum number of digits.

mem_p = digits_p
mem_q = digits_q
mem_o = digits_o

# A standard multiplication algorithm requires the inputs (p, q) and the output (o)
# to be stored in memory simultaneously. Therefore, the minimum total memory
# is the sum of the memory for each.
m = mem_p + mem_q + mem_o

# As requested, printing each number in the final equation.
print(f"{mem_p} + {mem_q} + {mem_o} = {m}")