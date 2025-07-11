# This script determines the smallest genus of a closed surface that can contain a given surface.

# Step 1: Define the properties of the initial surface, Sigma.
# The problem states Sigma has a genus of 10.
genus_Sigma = 10
# The problem states Sigma has a single boundary component.
num_boundaries_Sigma = 1
# The problem also specifies the boundary is 'unknotted'.

# Step 2: Formulate the method for creating the closed surface, Sigma'.
# To create a closed surface (one with no boundaries) from Sigma, we must attach
# a "cap" surface, C, to the boundary of Sigma.
# The genus of the resulting surface, Sigma', is given by the formula:
# genus(Sigma') = genus(Sigma) + genus(C)

# Step 3: Find the cap C with the minimum possible genus.
# To find the smallest genus for Sigma', we must choose the cap C with the
# smallest possible genus. The cap C must have one boundary component to match Sigma.
# The simplest surface with one boundary component is a disk. A disk has genus 0.
genus_C_min = 0

# The condition that the boundary is 'unknotted' is what makes it possible to use
# a disk of genus 0 as the cap while keeping the resulting surface embedded in R^3.

# Step 4: Calculate the minimum genus 'g' for the final surface Sigma'.
g = genus_Sigma + genus_C_min

# Print the step-by-step reasoning and the final calculation.
print("--- Problem Analysis ---")
print(f"The initial surface, Sigma, has a genus of {genus_Sigma}.")
print("It has a single boundary component, which is an unknotted circle.")
print("\n--- Construction of the Closed Surface, Sigma' ---")
print("To form a closed surface, we cap the boundary of Sigma with another surface, C.")
print("The genus of the resulting surface follows the rule: g(Sigma') = g(Sigma) + g(C).")
print("\n--- Minimizing the Genus ---")
print("To achieve the smallest possible genus for Sigma', we must choose a cap C with the minimum genus.")
print("The simplest surface with one boundary is a disk, which has genus 0.")
print("This is a valid choice because the unknotted boundary of Sigma can be filled by a disk in R^3.")
print(f"Therefore, the minimum genus for the cap is g(C)_min = {genus_C_min}.")
print("\n--- Final Calculation ---")
print("The smallest possible genus 'g' for Sigma' is:")
print(f"g = g(Sigma) + g(C)_min")
print(f"g = {genus_Sigma} + {genus_C_min}")
print(f"g = {g}")
