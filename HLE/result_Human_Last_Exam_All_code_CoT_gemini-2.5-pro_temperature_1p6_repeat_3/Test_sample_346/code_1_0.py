import sys
import math

# The evaluation of complex SKI expressions can lead to deep recursion.
# We increase Python's recursion limit to prevent any errors.
sys.setrecursionlimit(10000)

# Define the fundamental S, K, and I combinators as higher-order functions.
S = lambda x: lambda y: lambda z: x(z)(y(z))
K = lambda x: lambda y: x
I = lambda x: x

# For building the main expression, we need some common compound combinators.
# B is the composition combinator: B f g x = f(g(x)).
B = S(S(K(S))(K))

# U is the self-application combinator: U x = x(x).
U = S(I)(I)

# c2 is the Church numeral for 2: c2 f x = f(f(x)).
c2 = S(B)(I)

# Now, we construct the full expression from the problem statement.
# The structure is S(I)(S(I)(S(I)(A)(B))), where:
# A = K(c2)
# B = S(B)(U)
# Due to Python's left-associative evaluation, this will be handled as (...(A))(B).
# The K combinator in A will discard B, so the problematic U term is never used.
full_expression = S(I)(S(I)(S(I)(K(c2))(S(B)(U))))

# To find the integer value of a Church numeral, we apply it to a successor
# function (like incrementing by 1) and a starting value (like 0).
increment = lambda value: value + 1
zero = 0

# This evaluation will yield the integer n.
n = full_expression(increment)(zero)

# Calculate the base-2 logarithm of n.
log2_n = math.log2(n)

# Print the final result in the requested format.
# Remember to output each number in the final equation!
print(f"The integer n is: {n}")
print(f"The equation is: log_2({n}) = {log2_n}")
