import sys
import math

# The calculation involves deep recursion as it simulates Church numerals.
# A high recursion limit is necessary for a potentially large number n.
sys.setrecursionlimit(3000)

# 1. Define the primitive S, K, and I combinators.
def S(x):
    """The S combinator: S x y z = x z (y z)"""
    return lambda y: lambda z: x(z)(y(z))

def K(x):
    """The K combinator: K x y = x"""
    return lambda y: x

def I(x):
    """The I combinator: I x = x"""
    return x

# 2. Construct the full expression.
# The combinator for composition B = S(KS)K. It can be shown that the
# expression S(K(S))(K) used in the problem is also a valid B combinator.
B_combinator = S(K(S))(K)

# The combinator for the successor function is SUCC = S(B).
SUCC = S(B_combinator)

# The Church numeral for ZERO is K(I).
ZERO = K(I)

# The expression contains two main sub-expressions, let's call them A and B.
# A = SUCC(I), where I is the Church numeral for 1. So, A is the numeral for 2.
A = SUCC(I)

# The other sub-expression involves S(I)(I).
SII = S(I)(I)

# B = SUCC(S(I)(I))
B = SUCC(SII)

# Now, we build the full combinator for the target Church numeral n.
church_n_combinator = S(I)(S(I)(S(I)(K(A)))(B))

# 3. Evaluate the combinator expression to get the specific Church numeral for n.
# We do this by applying it to the SUCC and ZERO combinators.
n_combinator = church_n_combinator(SUCC)(ZERO)

# 4. Convert the resulting combinator n_combinator into a regular integer.
# We apply it to a simple increment function and the integer 0.
n = n_combinator(lambda i: i + 1)(0)

# 5. Calculate log_2 of the resulting integer n.
log_n = math.log2(n)

# 6. Print the results in the requested format.
print("The expression represents a Church integer n, which we can calculate.")
print(f"The integer value n is: {n}")
print("The problem asks for log_2(n).")
print(f"The equation is: log_2({n}) = {log_n}")
