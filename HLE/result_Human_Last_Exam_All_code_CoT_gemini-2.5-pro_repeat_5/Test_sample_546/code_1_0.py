# The value of p provided in the problem
p = 18446744074401676349

# As derived from the number theory analysis, the entire complex expression
# simplifies to 2^81, regardless of the specific large value of p,
# as long as the primality conditions hold.
base = 2
exponent = 81

# Calculate the final result
result = base ** exponent

# Print the final equation with all the numbers
print(f"The value of the function f(p) is calculated as follows:")
print(f"f({p}) = {base}^{exponent}")
print(f"f({p}) = {result}")
