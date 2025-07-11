import math
import sys

# Set a higher recursion limit for safety, as combinator logic can involve deep recursion.
sys.setrecursionlimit(2000)

# Define the basic SKI combinators as Python lambda functions.
S = lambda x: lambda y: lambda z: x(z)(y(z))
K = lambda x: lambda y: x
I = lambda x: x

# The given SKI combinator expression represents a Church numeral n.
# We build it by directly translating the expression into nested function calls in Python.
# Expression: S(I)(S(I)(S(I)(K(S(S(K(S))(K))(I)))))(S(S(K(S))(K))(S(I)(I))))

church_numeral_n = S(I)(
    S(I)(
        S(I)(
            K(
                S(S(K(S))(K))(I)
            )
            (
                S(S(K(S))(K))(S(I)(I))
            )
        )
    )
)

# A Church numeral n is a function that applies a given function f, n times to a value x.
# To find the integer value of n, we can apply it to a successor function (lambda i: i + 1)
# and a starting value (0). The result will be n.
successor_function = lambda x: x + 1
starting_value = 0

n = church_numeral_n(successor_function)(starting_value)

# Now, we calculate log_2(n) as requested.
result = math.log2(n)

# Print the determined value of n and the final equation.
print(f"The SKI expression evaluates to the Church numeral n = {n}.")
print(f"The final equation is: log_2({int(n)}) = {result}")
