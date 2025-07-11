# Define the numerical components of the scaling equation derived from the epsilon expansion.
# This is for a scalar field theory (N=1 universality class) to first order in epsilon.
numerator = 1
denominator = 6
upper_critical_dimension = 4

# Print the final quantitative scaling relationship.
# The formula expresses the specific heat critical exponent 'α'
# as a function of the spatial dimensionality 'd' for d < 4.
print("In scalar field theory, the scaling of the specific heat exponent α with dimensionality d (for d < 4) is given by:")
print(f"α = ({numerator}/{denominator}) * ({upper_critical_dimension} - d)")
