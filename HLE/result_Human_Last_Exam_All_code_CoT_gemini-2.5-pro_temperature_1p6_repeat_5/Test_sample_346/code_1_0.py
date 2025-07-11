import math

# Define the fundamental SKI combinators as lambda functions.
# S = λx.λy.λz. xz(yz)
S = lambda x: lambda y: lambda z: x(z)(y(z))
# K = λx.λy. x
K = lambda x: lambda y: x
# I = λx. x
I = lambda x: x

# We can build more complex combinators from these.
# For example, the Church numeral for 2 is S(I)(I)
# Let's construct the full expression given in the problem statement.
# To avoid confusion, we will build it piece by piece from the inside out.
# Application is left-associative: f(x)(y) is ((f(x))(y)).

# The core operator: Op = S(S(K(S))(K))
op_s = S(S(K(S))(K))

# The two main arguments in the middle of the expression:
# Arg 1 is K(Op(I))
term4_arg = op_s(I)
term4 = K(term4_arg)

# Arg 2 is Op(S(I)(I))
church_2 = S(I)(I)
term5 = op_s(church_2)

# The successor function: SUCC = S(I)
succ = S(I)

# The expression has the form SUCC( SUCC( (SUCC(term4))(term5) ) )
# We build it step-by-step
sub_expr_1 = succ(term4)
sub_expr_2 = sub_expr_1(term5)
sub_expr_3 = succ(sub_expr_2)
final_expression = succ(sub_expr_3)

# The final_expression is a Church integer 'n'.
# A Church integer is a function that takes two arguments, a function `f` and a value `x`,
# and applies `f` to `x`, `n` times.
# We can find the integer value of `n` by providing Python's increment function and 0.
py_succ = lambda i: i + 1
py_zero = 0

# Calculate the integer value of n
n = final_expression(py_succ)(py_zero)

# Calculate the base-2 logarithm of n
result = math.log2(n)

print(f"The SKI expression evaluates to the Church integer n = {n}.")
print(f"The equation we want to solve is log_2(n).")
print(f"Substituting the value of n, we get:")
print(f"log_2({n}) = {result}")