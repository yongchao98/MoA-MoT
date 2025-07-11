import numpy as np

# Set a random seed for reproducibility
np.random.seed(0)

print("Demonstration: Why an arbitrary vector beam cannot be generated from a fixed linear polarization input.")
print("-" * 70)

# 1. Define the random medium at a point via a random 2x2 Jones Matrix T.
# This matrix describes how polarization is transformed at this point.
T = np.random.rand(2, 2) + 1j * np.random.rand(2, 2)
print(f"Random medium's Jones matrix T:\n{T}\n")

# 2. Define the input beam with a fixed horizontal linear polarization.
# The controllable part is a scalar amplitude/phase U_in. We'll use U_in=1 for simplicity.
U_in = 1.0
E_in_polarization = np.array([1, 0])
E_in = U_in * E_in_polarization
print(f"Input beam's fixed polarization vector: {E_in_polarization}\n")

# 3. Calculate the output field's polarization state after passing through T.
# The output polarization is given by the first column of T, scaled by U_in.
E_out_vector = T @ E_in
print(f"The resulting output polarization vector is proportional to:\n{E_out_vector}\n")

# 4. Define an arbitrary target polarization, for example, right-hand circular polarization.
E_target_polarization = np.array([1, 1j]) / np.sqrt(2)
print(f"An arbitrary target polarization vector (e.g., circular):\n{E_target_polarization}\n")

# 5. Analyze the constraint by comparing the ratio of the vector components.
# To generate the target, the ratio of y/x components must match.
ratio_out = E_out_vector[1] / E_out_vector[0]
ratio_target = E_target_polarization[1] / E_target_polarization[0]

print("--- Analysis ---")
print(f"The ratio of (y/x) components for the output is fixed by T: {ratio_out:.4f}")
print(f"The ratio of (y/x) components for the arbitrary target is: {ratio_target:.4f}")
print("\nSince we can only control a single scalar value (U_in), we cannot independently")
print("set the two components of the output vector to match an arbitrary target.")
print("The direction of the output polarization is fixed by the medium.")
print("\nTherefore, the answer to the question is NO.")
