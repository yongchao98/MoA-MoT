import math

# The question asks for the value of the "computational factor" or "mushy zone constant"
# used in the original simulation-only paper that introduced the enthalpy-porosity method
# for phase change simulation.

# This method was formally presented in the paper by V. R. Voller and C. Prakash in 1987.
# In their paper, they defined a momentum source term for the mushy zone:
# S = -C * ((1-lambda)^2 / (lambda^3 + b)) * u
# where C is the mushy zone constant.

# The value they chose for this constant in their work is the answer to the question.
# This value was later modified in a subsequent paper by Brent, Voller, and Reid (1988)
# for a specific benchmark against gallium melting experiments.

# From the original Voller and Prakash (1987) paper, the value is specified.
original_computational_factor = 1.6 * math.pow(10, 6)

# The final equation is simply the value of this constant.
# We will print the components of this value as requested.
base = 1.6
exponent = 6

print(f"The computational factor used in the prior published implementation of this numerical method was C = {base} x 10^{exponent}.")
print(f"This is equal to: {int(original_computational_factor)}")
