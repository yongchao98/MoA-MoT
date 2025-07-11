import math

# Step 1: Define the given parameters.
# The problem gives k1 and k2 in scientific notation, which can be very large.
# We will primarily work with their base-10 logarithms to manage the scale.
log10_k1 = 3000
log10_k2 = 500
# lambda_2 = 0.5 * lambda_1 = 0.5, which implies lambda_1 = 1.0 and lambda_2 = 0.5
lambda_1 = 1.0
lambda_2 = 0.5
# The norm of h
h_norm = 1000
log10_h_norm = math.log10(h_norm)

# Step 2: Calculate the value for limsup ||x_n|| as n -> +infinity.
# According to Theorem 3.2 from the referenced paper,
# limsup_{n->+inf} ||x_n|| <= k1 * k2 * (1 - lambda2)^(-1) * |||h|||.
# We assume this bound is attained.
# To compute this, we use logarithms to handle the large numbers.
log10_L_plus_inf = log10_k1 + log10_k2 + math.log10(1 / (1 - lambda_2)) + log10_h_norm

# The actual value of L_plus_inf would be 10**log10_L_plus_inf.
# L_plus_inf = (10**log10_k1) * (10**log10_k2) * (1 / (1 - lambda_2)) * h_norm
# L_plus_inf = 10**3000 * 10**500 * 2 * 1000 = 2 * 10**3503
# For printing purposes, we can represent it as a string.
L_plus_inf_str = f"2e{int(log10_k1 + log10_k2 + log10_h_norm)}"

# Step 3: Calculate the value for liminf ||x_n|| as n -> -infinity.
# According to Theorem 3.3 from the paper,
# limsup_{n->-inf} ||x_n|| <= k2 * (1 - lambda2)^(-1) * |||h|||.
# We assume the limit exists, so liminf = limsup, and that the bound is attained.
log10_L_minus_inf = log10_k2 + math.log10(1 / (1 - lambda_2)) + log10_h_norm

# The actual value of L_minus_inf would be 10**log10_L_minus_inf.
# L_minus_inf = (10**log10_k2) * (1 / (1 - lambda_2)) * h_norm
# L_minus_inf = 10**500 * 2 * 1000 = 2 * 10**503
# For printing purposes, we can represent it as a string.
L_minus_inf_str = f"2e{int(log10_k2 + log10_h_norm)}"


# Step 4: Calculate the two terms in the final expression.
# The expression is: 100 * limsup log10(1/3 * ||x_n||) + 10 * liminf log10(1/3 * ||x_n||)
# Since log is continuous, we can swap lim and log.
# First term: 100 * log10(1/3 * L_plus_inf)
term1 = 100 * (math.log10(1/3) + log10_L_plus_inf)

# Second term: 10 * log10(1/3 * L_minus_inf)
term2 = 10 * (math.log10(1/3) + log10_L_minus_inf)

# Step 5: Calculate the final result.
result = term1 + term2

# Step 6: Print the output in the requested format.
# We need to print the numbers in the final equation.
# The large numbers L_plus_inf and L_minus_inf are represented by strings.
print(f"The expression to calculate is:")
print(f"100 * log10(1/3 * ||x_n||) for n->+inf   +   10 * log10(1/3 * ||x_n||) for n->-inf")
print("\nAssuming the bounds from the reference are sharp and the limit at -inf exists, we have:")
print(f"limsup ||x_n|| as n->+inf = {L_plus_inf_str}")
print(f"liminf ||x_n|| as n->-inf = {L_minus_inf_str}")
print("\nSubstituting these values, the calculation is:")
print(f"100 * log10(1/3 * {L_plus_inf_str}) + 10 * log10(1/3 * {L_minus_inf_str}) = {result}")
