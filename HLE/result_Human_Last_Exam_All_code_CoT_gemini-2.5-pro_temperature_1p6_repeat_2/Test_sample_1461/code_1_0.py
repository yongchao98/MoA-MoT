# Plan: To determine the helical pattern for an alternating alpha/epsilon-amino acid foldamer,
# we will extrapolate from known trends of similar foldamers.

# Step 1: Define the amino acids based on the number of carbons (p) between the N and C termini.
# This value determines the length of the backbone.
p_gamma = 3
p_delta = 4
p_epsilon = 5

print(f"Let p be the number of backbone carbons separating the amino and carboxyl groups in the non-alpha monomer.")
print(f"For a gamma-amino acid, p = {p_gamma}")
print(f"For a delta-amino acid, p = {p_delta}")
print(f"For our epsilon-amino acid, p = {p_epsilon}\n")


# Step 2: Establish the known helical patterns for shorter alternating foldamers.
# The pattern is noted as x/y.
helix_alpha_gamma = (11, 13)
helix_alpha_delta = (12, 14)

print(f"The known helical pattern for alpha/gamma-peptides (p={p_gamma}) is {helix_alpha_gamma[0]}/{helix_alpha_gamma[1]}.")
print(f"The known helical pattern for alpha/delta-peptides (p={p_delta}) is {helix_alpha_delta[0]}/{helix_alpha_delta[1]}.\n")

# Step 3: Identify the linear trend.
# When p increases by 1 (from gamma to delta), the 'x' and 'y' values of the helix also increase by 1.
x_gamma = helix_alpha_gamma[0]

# Step 4: Extrapolate the trend to predict the pattern for the alpha/epsilon-peptide.
# We predict the 'x' value for the epsilon-foldamer based on the value for the gamma-foldamer.
p_diff = p_epsilon - p_gamma
predicted_x = x_gamma + p_diff

# The 'y' value in the x/y notation has consistently been x + 2.
predicted_y = predicted_x + 2

print("We can observe a linear trend. We will now extrapolate to predict the pattern for alpha/epsilon.")
print("The calculation for the first number in the pattern (x) is:")
print(f"x_epsilon = x_gamma + (p_epsilon - p_gamma)")
print(f"x_epsilon = {x_gamma} + ({p_epsilon} - {p_gamma})")
print(f"x_epsilon = {predicted_x}\n")

print("The calculation for the second number in the pattern (y) is:")
print(f"y_epsilon = x_epsilon + 2")
print(f"y_epsilon = {predicted_x} + {2}")
print(f"y_epsilon = {predicted_y}\n")

print(f"Therefore, the predicted helical pattern for an alpha/epsilon-amino acid foldamer is: {predicted_x}/{predicted_y}")
