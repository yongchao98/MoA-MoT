import math

def double_factorial(n):
    """Calculates the double factorial n!!."""
    if n < 0:
        raise ValueError("Double factorial is not defined for negative numbers.")
    if n == 0 or n == -1:
        # By convention, 0!! = 1 and (-1)!! = 1
        return 1
    result = 1
    for i in range(n, 0, -2):
        result *= i
    return result

# The sequence a(n) is for the number of non-vanishing Feynman diagrams
# of order 2n. We want to find a(3).
# The recurrence relation is:
# a(n) = (2n-1)!! - sum_{k=1..n-1} a(k) * (2(n-k)-1)!!
# The base case is a(0) = 1.

# Store computed values of a(n)
a = {}

# a(0) is the base case
a[0] = 1

# Calculate a(1)
# a(1) = (2*1-1)!! = 1!!
a[1] = double_factorial(1)

# Calculate a(2)
# a(2) = (2*2-1)!! - a(1) * (2*(2-1)-1)!!
# a(2) = 3!! - a(1) * 1!!
a[2] = double_factorial(3) - (a[1] * double_factorial(1))

# Now, calculate a(3) using the previously computed values
# a(3) = (2*3-1)!! - [a(1)*(2*(3-1)-1)!! + a(2)*(2*(3-2)-1)!!]
# a(3) = 5!! - [a(1)*3!! + a(2)*1!!]
df5 = double_factorial(5)
df3 = double_factorial(3)
df1 = double_factorial(1)
term1 = a[1] * df3
term2 = a[2] * df1
a[3] = df5 - (term1 + term2)

# Print the final equation with all numbers filled in.
print(f"To find a(3), we first need a(1) and a(2):")
print(f"a(1) = 1!! = {a[1]}")
print(f"a(2) = 3!! - (a(1) * 1!!) = {df3} - ({a[1]} * {df1}) = {a[2]}")
print("\nNow we can calculate a(3):")
print(f"a(3) = 5!! - (a(1) * 3!! + a(2) * 1!!)")
print(f"a(3) = {df5} - ({a[1]} * {df3} + {a[2]} * {df1})")
print(f"a(3) = {df5} - ({term1} + {term2})")
print(f"a(3) = {df5} - {term1 + term2}")
print(f"a(3) = {a[3]}")
