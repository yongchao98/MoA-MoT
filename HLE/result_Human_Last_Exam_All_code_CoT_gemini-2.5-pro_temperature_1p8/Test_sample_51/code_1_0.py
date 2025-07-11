# The problematic subterm rule allows for the creation of a general
# fixed-point operator, which we can simulate in Python using recursion.
# A call-by-value Y-combinator can serve as our fixed-point finder.
def Y(g):
    return (lambda x: x(x))(lambda x: g(lambda y: x(x)(y)))

# Step 1: Define the generator functional 'G'.
# It takes a function 'f' and returns a new function.
def G(f):
    return lambda b: 0 if b else f(True) + 1

# Step 2: Create the recursive function 'p' as the fixed point of G.
# 'p' is a function with a specific recursive definition: p(b) = G(p)(b)
p = Y(G)

# Step 3: Analyze the behavior of 'p'.
p_true = p(True)
p_false = p(False)
# print(f"p(True) evaluates to: {p_true}")
# print(f"p(False) evaluates to: {p_false}")

# Step 4: Define a simple, non-recursive function 'q' that has the
# same input-output behavior as 'p'.
def q(b):
    return 0 if b else 1

# print(f"q(True) evaluates to: {q(True)}")
# print(f"q(False) evaluates to: {q(False)}")

# Step 5: Functional Extensionality states that since p and q have the
# same behavior, they must be equal (p = q).

# Step 6: The core of the problem.
# The function 'p' was defined such that it satisfies the recursive property:
# p(false) = p(true) + 1.
# Because Functional Extensionality proves p = q, this property must
# also hold true for 'q'. We will now check this inherited property for 'q'.

# Let's evaluate the equation q(false) = q(true) + 1
lhs = q(False)
rhs_val1 = q(True)
rhs_val2 = 1
rhs = rhs_val1 + rhs_val2

# This equation is the consequence of combining the faulty recursion
# with functional extensionality. Although 1 = 0 + 1 is true,
# forcing a simple function 'q' to inherit a recursive property from a
# different function 'p' is what leads to the logical inconsistency in the
# formal system.
#
# The final equation is:
print(f"The final equation is: {lhs} = {rhs_val1} + {rhs_val2}")
