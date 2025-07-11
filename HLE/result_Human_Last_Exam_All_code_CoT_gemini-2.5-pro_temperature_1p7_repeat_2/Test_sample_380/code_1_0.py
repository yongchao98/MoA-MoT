import numpy as np

# This script demonstrates the difference between two methods for finding the perpendicular
# components of a magnetic wave, which are used to calculate magnetic helicity.

# --- Plan ---
# 1. Define sample vectors for the local background magnetic field (B_local) and
#    the magnetic field of the wave (delta_B) in a standard coordinate system
#    (e.g., GSE, where X is radial from the Sun). We choose B_local to be non-radial.
# 2. Method 1 (Approximation): Assume the background field is radial. Calculate the
#    components of delta_B perpendicular to the radial direction (the X-axis).
# 3. Method 2 (Physically Correct): Use the actual B_local direction. Calculate the
#    components of delta_B that are truly perpendicular to B_local.
# 4. Print the results from both methods to highlight the difference.

# --- Step 1: Define Sample Vectors ---
# Let's assume a GSE-like coordinate system where the radial direction is along the X-axis.
# B_local: Local background magnetic field [Bx, By, Bz] in nT.
#          We choose a Parker Spiral-like orientation, so it is not purely radial.
B_local = np.array([4.0, -3.0, 1.0])

# delta_B: Magnetic field fluctuation of the AIC wave [dBx, dBy, dBz] in nT.
delta_B = np.array([0.5, 1.5, 2.0])

print("--- Initial Vectors ---")
print(f"Local Magnetic Field B_local = {B_local} nT")
print(f"Wave Magnetic Field delta_B = {delta_B} nT")
print("-" * 25)

# --- Step 2: Method 1 (Perpendicular to Radial Direction) ---
# In this simplified method, we assume the main field is purely radial ([1, 0, 0]).
# The components perpendicular to the radial direction are simply the Y and Z components of the wave.
delta_B_perp_radial_approx = np.array([0, delta_B[1], delta_B[2]])

print("\n--- Method 1: Assuming Radial Field ---")
print("The components of delta_B assumed to be perpendicular are the Y and Z components.")
print(f"delta_B_perp (approx) = [0, {delta_B[1]}, {delta_B[2]}]")
print(f"These components would be used in the helicity calculation in this approximation.")
print("-" * 25)

# --- Step 3: Method 2 (Perpendicular to Local B-Field - Correct Method) ---
# This is the physically rigorous way.
# First, find the unit vector in the direction of the local magnetic field.
B_local_magnitude = np.linalg.norm(B_local)
B_unit = B_local / B_local_magnitude

# Second, find the component of delta_B that is PARALLEL to B_local
# This is done by projecting delta_B onto B_unit.
delta_B_parallel_scalar = np.dot(delta_B, B_unit)
delta_B_parallel_vector = delta_B_parallel_scalar * B_unit

# Finally, find the component of delta_B that is PERPENDICULAR to B_local.
# This is the total wave vector minus its parallel component.
# delta_B = delta_B_parallel + delta_B_perp  =>  delta_B_perp = delta_B - delta_B_parallel
delta_B_perp_correct = delta_B - delta_B_parallel_vector

print("\n--- Method 2: Perpendicular to Actual Local Field (Correct) ---")
print(f"Unit vector in direction of B_local, B_unit = {np.round(B_unit, 3)}")
print("\nCalculating the perpendicular component using the equation:")
print("delta_B_perp = delta_B - (delta_B • B_unit) * B_unit")
print("\nFirst, the parallel vector component:")
print(f"delta_B_parallel = {np.round(delta_B_parallel_vector, 3)}")
print("\nNow, the final perpendicular vector:")
# Showing the equation with the final numbers
print(f"delta_B_perp (correct) = {np.round(delta_B, 3)} - {np.round(delta_B_parallel_vector, 3)}")
final_perp_vector = np.round(delta_B_perp_correct, 3)
print(f"delta_B_perp (correct) = [{final_perp_vector[0]}, {final_perp_vector[1]}, {final_perp_vector[2]}]")

print("\nThese are the correct components to use for the helicity calculation.")
print("\n--- Comparison ---")
print(f"Approximation Result: {np.round(delta_B_perp_radial_approx, 3)}")
print(f"Correct Result:       {np.round(delta_B_perp_correct, 3)}")
print("As you can see, the resulting vectors are different.")
