# Define string representations for all variables in the model
# to construct the final equation.
b = "b"
c = "c"
pg = "pg"
pt = "pt"
gamma_t = "𝛾t"
mu_t = "𝜇t"
mu_g = "𝜇g"

# The expression for R0f is the product of two components:
# 1. The number of grass patches an average burning tree ignites.
# 2. The number of trees an average burning grass patch ignites.
# The final equation combines these two parts.

print("The expression for R0f is:")
# The f-string below constructs and prints the final equation,
# showing how each variable contributes.
print(f"R0f = ({b} * {c} * {pg} * {pt}) / (({gamma_t} + {mu_t}) * {mu_g})")
