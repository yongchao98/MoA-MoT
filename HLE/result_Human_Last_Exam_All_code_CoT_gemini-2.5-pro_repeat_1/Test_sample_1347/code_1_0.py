# Define string representations for the variables to use in the output.
R0f = "R0f"
b = "b"
pg = "pg"
gamma_t = "𝛾t"
mu_t = "𝜇t"
tau = "𝜏"
mu_g = "𝜇g"
c = "c"
pt = "pt"

# The final expression is constructed from several parts.
# We will define and print each part to clearly show how the final equation is built.

# Part 1: Number of grass areas ignited by one tree over its lifetime
grass_ignited_numerator = f"{b} * {pg}"
grass_ignited_denominator = f"({gamma_t} + {mu_t})"
grass_ignited_term = f"({grass_ignited_numerator}) / ({grass_ignited_denominator})"

# Part 2: Expected number of trees ignited by one grass area
prob_grass_survives_latency = f"{tau} / ({tau} + {mu_g})"
trees_ignited_by_infectious_grass = f"({c} * {pt}) / {mu_g}"
trees_ignited_term = f"({prob_grass_survives_latency}) * ({trees_ignited_by_infectious_grass})"

# Combine the parts for the final expression
final_numerator = f"{b} * {pg} * {c} * {pt} * {tau}"
final_denominator = f"({gamma_t} + {mu_t}) * ({tau} + {mu_g}) * {mu_g}"

print("The expression for R0f is the product of two terms:")
print("1. The number of grass areas a single burning tree ignites during its lifetime.")
print("2. The expected number of trees a single ignited grass area ignites during its lifetime.")
print("-" * 50)

print(f"Term 1 (Grass ignitions per tree): {grass_ignited_term}")
print(f"Term 2 (Tree ignitions per grass area): {trees_ignited_term}")
print("-" * 50)

print("By multiplying these terms, we get the final expression for R0f.")
print("The components of the final equation are broken down as follows:")
print(f"  Numerator (factors promoting fire spread): {final_numerator}")
print(f"  Denominator (factors limiting fire spread): {final_denominator}")
print("-" * 50)

print("The final expression for R0f is:")
print(f"{R0f} = ({final_numerator}) / ({final_denominator})")