# The problem asks for the value of a computational factor used in the original
# simulation-only implementation of the enthalpy-porosity method for melting.

# This method was detailed in the paper:
# Brent, A. D., Voller, V. R., and Reid, K. J. (1988).
# "Enthalpy-porosity technique for modeling convection-diffusion phase change:
# Application to the melting of a pure metal."
# Numerical Heat Transfer, Part A: Applications, 13(3), 297-318.

# In this paper, the momentum source term is based on the Carman-Kozeny equation
# and is controlled by a "mushy zone constant," denoted as C. This is the
# "computational factor" referred to in the problem.

# The paper explicitly states the value used for this constant in their simulations.
# We will define the base and exponent for this value.

base = 1.6
exponent = 6

# The value of the constant C is base * 10^exponent.
value = base * (10**exponent)

# Now, we print the final equation showing each number.
print(f"The value of the computational factor (C) from the prior published work is:")
print(f"C = {base} * 10^{exponent}")

# The numerical result is 1,600,000.
# This corresponds to option B.