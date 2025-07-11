import sympy as sp

# Define the symbols for the functions and variables
t = sp.Function('t')
F = sp.Function('F')
p = sp.Function('p')
theta = sp.Symbol('theta')

# The types of firms
theta_malpractice = 0
theta_truthful = 1

# The types of auditors
x_lenient = 0
x_strict = 1

# The derived condition for the sorting where more truthful firms (theta=1)
# choose more lenient auditors (x=0) and less truthful firms (theta=0)
# choose stricter auditors (x=1) is C(0, 1) < C(0, 0).
# C(0, 1) = p(1) * F(0)
# C(0, 0) = t(0)
# So the condition is t(0) > p(1) * F(0).

# We will print this condition.
condition_lhs = t(theta_malpractice)
condition_rhs = p(x_strict) * F(theta_malpractice)

# The instruction is to output each number in the final equation.
# Here we will output the symbolic representation.
print(f"The required condition is that the tax liability for a malpractice firm, t(0), must be greater than the expected penalty from using a strict auditor.")
print("This can be written as:")
print(f"{condition_lhs} > {condition_rhs}")