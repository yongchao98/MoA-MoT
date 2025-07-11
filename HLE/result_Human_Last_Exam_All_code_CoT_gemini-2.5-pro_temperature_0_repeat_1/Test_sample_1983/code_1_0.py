import math

# Step 1: Define the given constants
k1 = 10**3000
k2 = 10**500
lambda1 = 0.5
lambda2 = 0.5
h_norm = 1000

# Step 2: Calculate the limit superior and limit inferior based on the interpretation
# L_plus_inf corresponds to lim sup ||x_n|| as n -> +inf
L_plus_inf = h_norm * k1 / (1 - lambda1)

# L_minus_inf corresponds to lim inf ||x_n|| as n -> -inf
# Based on our assumption, this equals lim sup ||x_n|| as n -> -inf
L_minus_inf = h_norm * k2 / (1 - lambda2)

# Step 3: Calculate the two terms of the expression
# First term: 100 * log10( (1/3) * ||x_n||_sup(+) )
arg1 = (1/3) * L_plus_inf
log1 = math.log10(1/3) + math.log10(2) + math.log10(h_norm) + 3000 # log10(L_plus_inf) = log10(2*1000*10^3000)
term1 = 100 * log1

# Second term: 10 * log10( (1/3) * ||x_n||_inf(-) )
arg2 = (1/3) * L_minus_inf
log2 = math.log10(1/3) + math.log10(2) + math.log10(h_norm) + 500 # log10(L_minus_inf) = log10(2*1000*10^500)
term2 = 10 * log2

# Step 4: Calculate the final result
result = term1 + term2

# Step 5: Print the breakdown of the calculation as requested
print("Calculation Breakdown:")
print("-" * 25)
print(f"Value for lim sup ||x_n|| as n->+inf (L+): {h_norm} * ({k1:.1e}) / (1 - {lambda1}) = {L_plus_inf:.2e}")
print(f"Value for lim inf ||x_n|| as n->-inf (L-): {h_norm} * ({k2:.1e}) / (1 - {lambda2}) = {L_minus_inf:.2e}")
print("-" * 25)
print("Final Equation: 100 * log10( (1/3) * L+ ) + 10 * log10( (1/3) * L- )")
print(f"= 100 * log10( (1/3) * {L_plus_inf:.2e} ) + 10 * log10( (1/3) * {L_minus_inf:.2e} )")
print(f"= 100 * log10( {arg1:.2e} ) + 10 * log10( {arg2:.2e} )")
print(f"= 100 * ({log1:.4f}) + 10 * ({log2:.4f})")
print(f"= {term1:.4f} + {term2:.4f}")
print(f"= {result:.4f}")
print("-" * 25)
print(f"Final Answer: {result}")
