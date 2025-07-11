# Step 1 & 3: Define fundamental rates based on our assumptions.
# To get a numeric answer, we can set one rate to a nominal value, like 1.
# From Assumption 1 (steady-state): rate of true speciation (lambda_c) = rate of true extinction (mu)
mu = 1.0
lambda_c = 1.0

# From Assumption 2 (symmetric pseudo-extinction): rate of anagenesis (lambda_a) = rate of pseudo-extinction from bifurcation
lambda_a = 0.5 * lambda_c

# Step 2: Formulate the extinction rates.
# Extinction rate for an Evolutionary Species (ES) is just the true extinction rate.
mu_E = mu

# Extinction rate for a Morphospecies (MS) is the sum of true extinction, 
# pseudo-extinction from anagenesis, and pseudo-extinction from bifurcating speciation.
mu_M = mu + lambda_a + 0.5 * lambda_c

# Step 4: Calculate the multiplicative factor.
factor = mu_M / mu_E

# Print the breakdown of the calculation and the final answer.
# The equation shows the components of the morphospecies extinction rate divided by the evolutionary species extinction rate.
print("The calculation of the multiplicative factor is:")
print(f"({mu} [true extinction] + {lambda_a} [anagenesis] + {0.5 * lambda_c} [bifurcation]) / {mu_E} [true extinction] = {factor}")
print("\nThe extinction rate for a morphospecies is greater than the extinction rate for an evolutionary species by a factor of:")
print(int(factor))
<<<2>>>