# The problem requires finding the original value of the "computational factor" 
# from the prior published simulation-only work on the enthalpy-porosity method.
# This factor is known as the mushy zone constant.
# The seminal "simulation-only" paper is Voller & Prakash (1987), "A fixed grid 
# numerical modelling methodology for convection-diffusion mushy region 
# phase-change problems."
# In that paper, the authors specify the value used for this constant in their work.

# Define the components of the value as per the scientific literature
base_value = 1.6
exponent = 6

# Calculate the final value
original_computational_factor = base_value * (10**exponent)

# The question asks to output the numbers in the final equation.
# Here, we reconstruct the value from its components to show the "equation".
print(f"The original computational factor is calculated as: {base_value} * 10^{exponent}")
print(f"Final Value: {original_computational_factor}")
