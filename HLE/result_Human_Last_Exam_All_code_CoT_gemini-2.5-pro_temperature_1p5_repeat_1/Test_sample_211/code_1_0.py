import numpy as np
from scipy.stats import linregress

def solve_diversity_slope():
    """
    Simulates the effect of an invasive species on the latitudinal diversity gradient
    under Hubbell's neutral theory.
    """
    # 1. Define the sites and their approximate latitudes
    sites = {
        "Tena, Ecuador": -0.99,
        "Cali, Colombia": 3.45,
        "Barro Colorado Island, Panama": 9.15,
        "Sarapiqui, Costa Rica": 10.46,
        "Managua, Nicaragua": 12.13,
        "Yoro, Honduras": 15.14,
        "Belmopan, Belize": 17.25,
        "Merida, Mexico": 20.97,
        "Miami, Florida": 25.76,
        "Charlotte, North Carolina": 35.23,
        "State college, Pennsylvania": 40.80,
        "Ottawa, Canada": 45.42
    }
    
    site_names = list(sites.keys())
    # Use absolute latitude for gradient calculation
    latitudes = np.array([abs(lat) for lat in sites.values()])

    # 2. Create a hypothetical initial diversity, showing a latitudinal gradient.
    # Model: Diversity decreases as we move away from the equator (latitude 0).
    # Let's say D = 250 - 4 * |latitude|
    initial_diversity = 250 - 4 * latitudes
    
    # 3. Calculate and print the initial slope
    slope_initial, intercept_initial, _, _, _ = linregress(latitudes, initial_diversity)
    
    print("--- Initial State: Pre-Invasion ---")
    print("A clear latitudinal diversity gradient exists.\n")
    for i in range(len(site_names)):
        print(f"{site_names[i]:<32}: Latitude {latitudes[i]:.2f}, Initial Diversity {initial_diversity[i]:.0f}")
    
    print("\nInitial slope of diversity vs. latitude:")
    print(f"y = {slope_initial:.2f}x + {intercept_initial:.2f}\n")
    
    # 4. Simulate biotic homogenization due to the invasive species.
    # Diversity at all sites collapses towards a low value (e.g., 5).
    # This models the invasive species out-competing most native species everywhere.
    final_diversity = (initial_diversity * 0.1) + (5 * 0.9) # Weighted average pulling all values towards 5
    
    # 5. Calculate and print the final slope
    slope_final, intercept_final, _, _, _ = linregress(latitudes, final_diversity)
    
    print("--- Final State: Post-Invasion ---")
    print("Biotic homogenization reduces and flattens diversity across all sites.\n")
    for i in range(len(site_names)):
        print(f"{site_names[i]:<32}: Latitude {latitudes[i]:.2f}, Final Diversity {final_diversity[i]:.0f}")
        
    print("\nFinal slope of diversity vs. latitude:")
    print(f"y = {slope_final:.2f}x + {intercept_final:.2f}\n")
    
    # 6. Conclusion
    print("--- Conclusion ---")
    print(f"The initial slope of {slope_initial:.2f} represents a strong diversity gradient.")
    print(f"After the invasion leads to homogenization, the slope flattens to {slope_final:.2f}, which is much closer to zero.")
    print("Ultimately, under Hubbell's unified theory, the slope of insect diversity will flatten and approach zero.")

solve_diversity_slope()