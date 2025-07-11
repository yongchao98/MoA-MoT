import sys
import math

# To handle the deep recursion that can occur when evaluating combinator expressions,
# we increase Python's default recursion limit to a sufficiently large value.
sys.setrecursionlimit(300000)

# Step 1: Define the fundamental SKI combinators as higher-order functions (lambdas).
S = lambda x: lambda y: lambda z: x(z)(y(z))
K = lambda x: lambda y: x
I = lambda x: x

# Step 2: Construct the complex SKI expression given in the problem.
# The expression is: S(I)(S(I)(S(I)(K(S(S(K(S))(K))(I)))))(S(S(K(S))(K))(S(I)(I))))
# We build this function by applying the combinators to each other, following the rules of application.

# Let's define the recurring sub-expression for clarity: A = S(S(K(S))(K))
A = S(S(K(S))(K))

# Another recurring part is the Church numeral for 2: C2 = S(I)(I)
C2 = S(I)(I)

# Now, we build the main expression from the inside out.
# The expression's structure is: S(I)(Y)(Z)
# where Y = S(I)(S(I)(K(A(I))))
# and   Z = A(C2)

# Build the Y part:
part_A_I = A(I)
part_K_A_I = K(part_A_I)
part_S_I_K_A_I = S(I)(part_K_A_I)
Y = S(I)(part_S_I_K_A_I)

# Build the Z part:
Z = A(C2)

# Assemble the final expression, which corresponds to our Church numeral 'n'.
church_n_func = S(I)(Y)(Z)

# Step 3: Evaluate the Church numeral to get its integer value.
# A Church numeral 'n' represents the function `lambda f: lambda x: f(f(...f(x)...))` (f applied n times).
# To find the integer value, we provide a successor function for 'f' and zero for 'x'.
successor = lambda i: i + 1
zero = 0

# The application church_n_func(successor)(zero) will compute the integer n.
n = church_n_func(successor)(zero)

# Step 4: Calculate log_2(n).
result = math.log2(n)

# Step 5: Print the final answer in the requested format, showing the calculated values.
print(f"log_2({n}) = {result}")
