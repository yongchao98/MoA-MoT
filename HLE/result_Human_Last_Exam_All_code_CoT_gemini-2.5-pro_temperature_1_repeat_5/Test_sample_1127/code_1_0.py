# The minimal polynomial for the connective constant (μ) is P(t) = t^4 - 4t^2 + 2.
# The equation is P(μ) = 0.
# The code below prints this equation, showing each non-zero coefficient as requested.

coeff_t4 = 1
coeff_t2 = -4
constant = 2

# We represent the equation with the variable 'μ' for the connective constant.
print("The final equation is:")
print(f"{coeff_t4} * μ^4 + ({coeff_t2}) * μ^2 + {constant} = 0")