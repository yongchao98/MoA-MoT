import numpy as np

def model_diversity_gradient():
    """
    This function models the change in the latitudinal diversity gradient
    of insects following an invasion, based on Hubbell's unified neutral theory.
    """
    # 1. Define the 12 study sites and their approximate latitudes (in degrees North)
    sites = {
        "Tena, Ecuador": 0.9,
        "Cali, Colombia": 3.4,
        "Barro Colorado Island, Panama": 9.1,
        "Sarapiqui, Costa Rica": 10.4,
        "Managua, Nicaragua": 12.1,
        "Yoro, Honduras": 15.0,
        "Belmopan, Belize": 17.2,
        "Merida, Mexico": 20.9,
        "Miami, Florida": 25.7,
        "Charlotte, North Carolina": 35.2,
        "State college, Pennsylvania": 40.8,
        "Ottawa, Canada": 45.4
    }

    latitudes = np.array(list(sites.values()))
    site_names = list(sites.keys())

    # 2. Simulate initial alpha diversity based on a latitudinal gradient.
    # We'll use a linear model: Diversity = intercept - slope * latitude.
    # Higher diversity at the equator (latitude ~0).
    initial_intercept = 220
    initial_slope = -3.5
    # Calculate initial diversity with some random noise for realism
    initial_diversity = initial_intercept + initial_slope * latitudes + np.random.normal(0, 5, len(latitudes))
    initial_diversity = initial_diversity.astype(int)

    print("--- Initial State: Pre-Invasion ---")
    print("A clear latitudinal diversity gradient exists.\n")
    
    # Fit a line to the initial data to find the slope
    fit_initial_slope, fit_initial_intercept = np.polyfit(latitudes, initial_diversity, 1)

    print(f"Initial Equation of Diversity vs. Latitude:")
    print(f"Diversity = {fit_initial_slope:.2f} * Latitude + {fit_initial_intercept:.2f}\n")

    # 3. Simulate the invasion effect according to Hubbell's theory.
    # The invasive species causes biotic homogenization, reducing diversity everywhere
    # and thus flattening the gradient.
    # We model this by reducing diversity across all sites, which lessens the
    # difference between tropical and temperate sites.
    
    # A simple way to model homogenization is to make all sites lose a
    # proportion of their species, which makes the absolute difference smaller.
    final_diversity = initial_diversity * 0.45 + np.random.normal(0, 3, len(latitudes))
    final_diversity = final_diversity.astype(int)

    print("--- Final State: Post-Invasion ---")
    print("The invasive species has reduced diversity, causing biotic homogenization.\n")

    # 4. Fit a line to the final data to find the new, flatter slope.
    fit_final_slope, fit_final_intercept = np.polyfit(latitudes, final_diversity, 1)

    print("Final Equation of Diversity vs. Latitude:")
    # We explicitly print each number for the final equation as requested
    slope_val = f"{fit_final_slope:.2f}"
    intercept_val = f"{fit_final_intercept:.2f}"
    print(f"Diversity = {slope_val} * Latitude + {intercept_val}\n")
    
    print("--- Conclusion ---")
    print(f"The initial slope was {fit_initial_slope:.2f}.")
    print(f"After the invasion, the new slope is {fit_final_slope:.2f}.")
    print("The slope has become less steep (flatter), moving closer to zero.")

# Run the simulation
model_diversity_gradient()