import numpy as np
from scipy import stats

def analyze_diversity_slope():
    """
    Analyzes the change in the slope of the latitudinal diversity gradient
    according to Hubbell's unified neutral theory after a widespread invasion.
    """

    # The 12 study sites and their approximate latitudes (in degrees)
    # Negative values represent the Southern Hemisphere.
    sites = {
        "Tena, Ecuador": -0.99,
        "Cali, Colombia": 3.45,
        "Barro Colorado Island, Panama": 9.15,
        "Sarapiqui, Costa Rica": 10.4,
        "Managua, Nicaragua": 12.13,
        "Yoro, Honduras": 15.14,
        "Belmopan, Belize": 17.25,
        "Merida, Mexico": 20.97,
        "Miami, Florida": 25.76,
        "Charlotte, North Carolina": 35.23,
        "State college, Pennsylvania": 40.80,
        "Ottawa, Canada": 45.42
    }

    latitudes = np.array(list(sites.values()))

    # --- Initial State: A clear Latitudinal Diversity Gradient ---
    # We will simulate a strong gradient where diversity is high at the equator
    # and low at the poles. This is just a hypothetical model for demonstration.
    # For example, let's say diversity drops by 4 species for every degree of latitude.
    initial_diversity = 180 - 3.5 * latitudes
    
    # Calculate the initial slope
    initial_slope, initial_intercept, _, _, _ = stats.linregress(latitudes, initial_diversity)

    print("--- Initial State (Before Widespread Invasion) ---")
    print(f"The initial slope of the diversity gradient is approximately: {initial_slope:.4f}")
    print("This steep negative slope indicates that species diversity is highest near the equator and decreases towards the poles.\n")

    # --- Final State: After the Invasive Species Dominates ---
    # According to Hubbell's neutral theory, ecological drift will ultimately lead
    # to the fixation of the single most abundant species (the invasive).
    # Therefore, the alpha diversity at every single site will become 1.
    final_diversity = np.ones_like(latitudes)

    # Calculate the final slope
    final_slope, final_intercept, _, _, _ = stats.linregress(latitudes, final_diversity)
    
    print("--- Final State (Long-Term Outcome of Invasion) ---")
    print("Ultimately, the invasive species will be the only one left at all sites.")
    print("The alpha diversity at every site will be 1.")
    print(f"The final slope of the diversity gradient is: {final_slope:.4f}")

    # As requested, printing each number in the final equation
    print("\nThe final equation describing the relationship between diversity and latitude will be:")
    print(f"Final Diversity = ({final_slope:.1f}) * Latitude + ({final_intercept:.1f})")


if __name__ == '__main__':
    analyze_diversity_slope()
