# The derived formula for f_n(p) contains several constant numbers.
# Let's define them to fulfill the instruction "output each number in the final equation".

# In the term (1-p)^n, the number is 1.
c1 = 1
# In the denominator (2p-1), the numbers are 2 and 1.
c2 = 2
c3 = 1

# Print the final formula using an f-string to place the numbers in context.
print("The simple formula for f_n(p) is:")
print(f"f_n(p) = (p^n - ({c1} - p)^n) / ({c2}*p - {c3})")
