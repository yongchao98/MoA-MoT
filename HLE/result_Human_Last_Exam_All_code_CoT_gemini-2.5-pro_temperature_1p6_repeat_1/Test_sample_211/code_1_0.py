import numpy as np

def calculate_and_display_slope(sites, diversity_values, state_name):
    """
    Calculates and displays the slope of a diversity gradient.
    The function uses numpy's polyfit to perform a linear regression.
    """
    # We use an index for the sites as our x-axis (from south to north)
    site_indices = np.arange(len(sites))
    
    # Use polyfit to find the slope (m) and intercept (b) of the line y = mx + b
    # The degree of the polynomial is 1 for a straight line
    slope, intercept = np.polyfit(site_indices, diversity_values, 1)
    
    print(f"--- {state_name} State ---")
    print(f"Site Names: {sites}")
    print(f"Diversity Values: {diversity_values}")
    print(f"The equation for the diversity gradient is approximately: Diversity = {slope:.2f} * (Site Index) + {intercept:.2f}")
    print(f"The slope of the diversity gradient is: {slope:.2f}\n")
    return slope, intercept

# --- 1. Define the sites and the initial state ---
# List of 12 sites from south (equator) to north
sites = [
    "Tena", "Cali", "Barro Colorado Island", "Sarapiqui",
    "Managua", "Yoro", "Belmopan", "Merida",
    "Miami", "Charlotte", "State college", "Ottawa"
]

# Hypothetical initial alpha diversity, showing a latitudinal gradient (high near equator, low near poles)
initial_diversity = [150, 145, 130, 115, 100, 90, 85, 70, 50, 40, 25, 15]

# --- 2. Define the final state, as predicted by neutral theory ---
# Ultimately, one species (the invader) dominates everywhere, making alpha diversity = 1 at all sites
final_diversity = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]


# --- 3. Calculate and display the results ---
# Calculate and show the slope for the initial state
calculate_and_display_slope(sites, initial_diversity, "Initial")

# Calculate and show the slope for the final state
final_slope, final_intercept = calculate_and_display_slope(sites, final_diversity, "Final")

# As requested, output the numbers in the final equation explicitly
print("--- Final Equation Breakdown ---")
print(f"Final Diversity = (Slope) * (Site Position) + (Intercept)")
print(f"Final Diversity = ({final_slope:.2f}) * (Site Position) + ({final_intercept:.2f})")

<<<0>>>