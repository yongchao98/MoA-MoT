# Define the positions of the portal and the torch
portal_bin = 2025
torch_bin = 2024

# The starting position of the marble
start_bin = 0

# The probability of jumping from bin n to bin n+i is (1/3)^|i|.
# We can think of the probability of being captured by the portal or the torch
# as being proportional to the probability of jumping directly to them.

# The "influence" or probability of a direct jump to the portal
prob_to_portal = f"(1/3)^{abs(portal_bin - start_bin)}"
val_prob_to_portal_exp = abs(portal_bin - start_bin)

# The "influence" or probability of a direct jump to the torch
prob_to_torch = f"(1/3)^{abs(torch_bin - start_bin)}"
val_prob_to_torch_exp = abs(torch_bin - start_bin)

# The probability of escaping is the ratio of the influence of the portal
# to the sum of the influences of both the portal and the torch.
# P(escape) = prob_to_portal / (prob_to_portal + prob_to_torch)
# P(escape) = (1/3)^2025 / ((1/3)^2025 + (1/3)^2024)
# We can simplify this by multiplying the numerator and denominator by 3^2025:
# P(escape) = 1 / (1 + 3^(2025-2024))
# P(escape) = 1 / (1 + 3)

numerator = 1
denominator = 1 + 3**(val_prob_to_portal_exp - val_prob_to_torch_exp)
result = numerator / denominator

print(f"The position of the escape portal is {portal_bin}.")
print(f"The position of the torch is {torch_bin}.")
print(f"The starting position of the marble is {start_bin}.")
print("The probability of escaping is calculated as the ratio of the probability of a direct jump to the portal to the sum of the probabilities of a direct jump to either the portal or the torch.")
print(f"Let P_escape be the probability of escaping and P_melt be the probability of melting.")
print(f"P_escape is proportional to (1/3)^|{portal_bin} - {start_bin}| = (1/3)^{val_prob_to_portal_exp}.")
print(f"P_melt is proportional to (1/3)^|{torch_bin} - {start_bin}| = (1/3)^{val_prob_to_torch_exp}.")
print(f"The desired probability is P_escape / (P_escape + P_melt).")
print(f"This simplifies to 1 / (1 + 3^({val_prob_to_portal_exp} - {val_prob_to_torch_exp})).")
print(f"So the final equation is {numerator} / ({numerator} + 3^({val_prob_to_portal_exp - val_prob_to_torch_exp})).")
print(f"The final probability is {numerator}/{denominator}, which is {result}.")