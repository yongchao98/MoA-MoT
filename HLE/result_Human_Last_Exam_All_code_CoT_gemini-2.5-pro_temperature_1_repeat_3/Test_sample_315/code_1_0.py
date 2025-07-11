import math

# Plan:
# 1. Define the physical parameters of the experiment.
# 2. Use the known specific rotation of D-glucose at a reference wavelength (Sodium D-line, 589 nm)
#    to calibrate the Drude equation, which models how specific rotation changes with wavelength.
# 3. Define the wavelengths for different representative colors (Red, Green, Blue).
# 4. For each color, calculate its specific rotation using the Drude equation.
# 5. For each color, calculate the total rotation angle after passing through the 1m tube using Biot's Law.
# 6. Print the results, showing that each color is rotated by a different amount, which is the
#    physical basis for the spiraling rainbow effect.

# --- Step 1: Define physical parameters ---
# Path length (l): 1 m = 10 dm
l = 10  # in decimeters (dm)
# Concentration (c): Assume a reasonably concentrated solution, e.g., 250 g/L = 0.25 g/mL
c = 0.25  # in g/mL

print("This script calculates the rotation of polarized light for different colors in a D-glucose solution.")
print(f"Path length (l) = {l} dm")
print(f"Concentration (c) = {c} g/mL\n")

# --- Step 2: Calibrate Drude Equation ---
# Specific rotation of D-glucose for Sodium D-line ([α]_D)
specific_rotation_ref = 52.7  # degrees * mL / (dm * g)
# Wavelength of Sodium D-line (λ_ref)
lambda_ref = 589  # in nm
# Characteristic wavelength for sugars in Drude's equation (λ₀)
lambda_0 = 150  # in nm

# Drude's one-term equation: [α] = A / (λ^2 - λ₀^2)
# We find the constant A using our reference value.
# A = [α]_ref * (λ_ref^2 - λ₀^2)
A_constant = specific_rotation_ref * (lambda_ref**2 - lambda_0**2)

print("--- Calculating Rotation for Different Colors ---\n")

# --- Step 3: Define colors and their wavelengths ---
colors = {
    "Red": 650,    # nm
    "Green": 550,  # nm
    "Blue": 450    # nm
}

# --- Steps 4-6: Calculate and Print Results ---
for color, lam in colors.items():
    # Step 4: Calculate specific rotation [α] for this wavelength
    specific_rotation_lambda = A_constant / (lam**2 - lambda_0**2)

    # Step 5: Calculate total rotation α using Biot's Law: α = [α] * l * c
    total_rotation = specific_rotation_lambda * l * c

    print(f"For {color} light (λ = {lam} nm):")
    # Step 6: Print the final equation with numbers
    print(f"Total Rotation α = [α] * l * c")
    print(f"α = {specific_rotation_lambda:.2f} * {l} * {c} = {total_rotation:.2f} degrees")
    print("-" * 20)

print("\nThe calculation shows that blue light is rotated significantly more than red light.")
print("This difference in rotation for each color is what causes the light to separate and form a spiraling rainbow when viewed from the side.")
