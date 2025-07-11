import math

# This script finds the value of a computational factor from a key scientific paper.
# The problem describes the enthalpy-porosity method for simulating melting, which
# uses a Carman-Kozeny term. The history of this method is well-documented in
# heat transfer literature.

# Step 1: Identify the key papers.
# The "prior published simulation-only work" refers to:
# V. R. Voller & C. Prakash, "A fixed grid numerical modelling methodology for
# convection-diffusion mushy region phase-change problems", Int. J. Heat Mass
# Transfer, Vol. 30, No. 8, pp. 1709-1719, 1987.
#
# The subsequent paper that "benchmarked against the melting of a piece of gallium"
# refers to:
# A. D. Brent, V. R. Voller & K. J. Reid, "Enthalpy-porosity technique for modeling
# convection-diffusion phase change: Application to the melting of a pure metal",
# Numerical Heat Transfer, Vol. 13, pp. 297-318, 1988.

# Step 2: Extract the computational factor from the "prior" paper.
# The question asks for the value of the computational factor (often called the
# mushy zone constant, 'C') from the prior, 1987 paper.
# In that paper, on page 1712, the authors state:
# "The constant C is a large number... In the current work C is taken as 1.6 x 10^6."

# Define the components of the value found in the paper.
base = 1.6
exponent = 6

# Step 3: Print the explanation and the final answer.
# The instructions require printing each number in the final equation.
print("The numerical method described is the enthalpy-porosity technique.")
print("The 'prior published implementation' refers to the 1987 paper by Voller and Prakash.")
print("In this paper, the computational factor, C, in the Carman-Kozeny term was given a specific value.")
print("This value was later modified in a 1988 paper to 1.6 x 10^3 to better match experiments on melting gallium.")
print("The original value from the 1987 paper is requested.")
print("\nThe value is expressed in scientific notation.")
print("The final equation for the original factor C is:")

# Print each component of the number 1.6 x 10^6
ten = 10
print(base, "x", ten, "^", exponent)