import math

# Step 1: Define the constants based on the problem description.
# Rs is the Schwarzschild Radius for Pandora. Since Pandora is 2x the Sun's mass,
# its Rs is 2x the Sun's Rs (~3 km).
Rs_pandora_km = 6

# r is the distance from the center of Pandora to the Pioneer probe.
# r = Pandora's radius (20 km) + distance from surface (60 km).
r_km = 80

# The divisor from the Taylor approximation for sqrt(1-x) ≈ 1 - x/2.
divisor = 2
one = 1

# Step 2: Calculate the gravitational time dilation factor using the approximation.
# f ≈ 1 - (Rs / r) / 2
dilation_factor = one - (Rs_pandora_km / r_km) / divisor

# Step 3: Determine the memory usage (z) for the most memory-efficient C program.
# The Wuxing 'frac' type size is 2D (n) + 2D (d) + 2D (e) = 6D.
# The most efficient program reuses a single 'frac' variable for all calculations.
memory_usage_D = 6

# Step 4: Format the output to include the numbers from the final equation
# and present the answer as f:z, as requested.
# The final equation is: f = 1 - (6 / 80) / 2
print(f"Based on the equation f = {one} - ({Rs_pandora_km} / {r_km}) / {divisor}, the result is {dilation_factor:.4f}:{memory_usage_D}")