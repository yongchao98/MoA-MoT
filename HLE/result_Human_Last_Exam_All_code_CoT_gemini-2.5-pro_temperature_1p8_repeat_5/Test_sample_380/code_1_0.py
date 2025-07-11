import numpy as np

# Plan:
# 1. Define a realistic, non-radial background magnetic field (B0) at L1, typical of the Parker spiral.
#    The coordinate system is assumed to be one where the X-axis is the Sun-Earth line (radial direction).
# 2. Construct a wave magnetic field (b_wave) that is perfectly perpendicular to B0,
#    as would be the case for an Alfvenic fluctuation.
# 3. Calculate the true power of this perpendicular wave: P_true = |b_wave|^2.
# 4. Calculate the approximated perpendicular power using only the Y and Z components: P_approx = b_wave_y^2 + b_wave_z^2.
# 5. Compare the two values to show the error introduced by the approximation.

print("--- Justification for Wave Coordinate System Choice ---\n")

# Step 1: Define the background magnetic field (B0) in nanoteslas (nT).
# At L1 (~1 AU), the Parker spiral angle is about 45 degrees relative to the radial line.
# A negative Bx represents a field pointing away from the Sun along the Sun-Earth line.
# We choose components Bx=-4 and By=4 to create this 45-degree angle.
B0 = np.array([-4.0, 4.0, 1.0])
print(f"Assumed background magnetic field B0 = {B0} nT")
print("Note: This field is NOT purely radial because the Y-component is non-zero.\n")

# Step 2: Construct a wave magnetic field (b_wave) that is perpendicular to B0.
# An easy way to find a vector perpendicular to B0 is to take the cross product of B0 with another vector.
# Let's use the Z-axis vector [0, 0, 1]. The resulting vector b_wave will be perpendicular to B0.
# For B0 = [Bx, By, Bz], the cross product B0 x [0, 0, 1] is [By, -Bx, 0].
b_wave_direction = np.array([B0[1], -B0[0], 0.0])

# We'll give the wave an amplitude of 1 nT for simplicity.
# To do this, we normalize the direction vector and multiply by the amplitude.
amplitude = 1.0
b_wave = (b_wave_direction / np.linalg.norm(b_wave_direction)) * amplitude
print(f"Constructed wave magnetic field b_wave = [{b_wave[0]:.3f}, {b_wave[1]:.3f}, {b_wave[2]:.3f}] nT")

# We can check that b_wave is perpendicular to B0. Their dot product should be zero.
dot_product = np.dot(b_wave, B0)
print(f"Sanity Check: Dot product of b_wave and B0 is {dot_product:.3f} (should be 0).\n")

# Step 3: Calculate the true perpendicular power of the wave.
# Wave power is proportional to the square of its amplitude. Since we constructed b_wave to be
# perfectly perpendicular to B0, its total power IS the true perpendicular power.
P_true = np.sum(b_wave**2)
print("--- Method 1: Rigorous Calculation (Field-Aligned) ---")
print("The wave's true perpendicular power is the square of its total magnitude.")
print("Equation: P_true = b_x² + b_y² + b_z²")
print(f"Calculation: P_true = {b_wave[0]:.3f}² + {b_wave[1]:.3f}² + {b_wave[2]:.3f}² = {P_true:.3f} nT²\n")

# Step 4: Calculate the approximated power using only Y and Z components.
# This simplified method incorrectly assumes the background field B0 is purely radial (along the X-axis).
P_approx = b_wave[1]**2 + b_wave[2]**2
print("--- Method 2: Simplified Calculation (Radial-Aligned) ---")
print("The simplified power calculation uses only components perpendicular to the radial (X) direction.")
print("Equation: P_approx = b_y² + b_z²")
print(f"Calculation: P_approx = {b_wave[1]:.3f}² + {b_wave[2]:.3f}² = {P_approx:.3f} nT²\n")

# Step 5: Compare the results and calculate the error.
error_percentage = 100 * abs(P_true - P_approx) / P_true
print("--- Comparison and Error Analysis ---")
print("The simplified method introduces an error by ignoring the wave's x-component.")
print("The final equation for the percentage error is:")
print(f"Error % = |(P_true - P_approx) / P_true| * 100")
print(f"Error % = |({P_true:.3f} - {P_approx:.3f}) / {P_true:.3f}| * 100 = {error_percentage:.1f}%")
