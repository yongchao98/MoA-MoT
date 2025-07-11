import numpy as np

def solve_diversity_slope():
    """
    This script explains and calculates the ultimate change in the slope of a 
    latitudinal diversity gradient after the introduction of a highly competitive 
    invasive species, based on the principles of Hubbell's Unified Neutral Theory.
    """

    # 1. Define the study sites and their approximate absolute latitudes.
    # The absolute value of latitude is used because diversity generally decreases 
    # in both directions away from the equator (latitude 0).
    sites = {
        "Tena, Ecuador": 0.99,
        "Cali, Colombia": 3.45,
        "Barro Colorado Island, Panama": 9.15,
        "Sarapiqui, Costa Rica": 10.46,
        "Managua, Nicaragua": 12.13,
        "Yoro, Honduras": 15.14,
        "Belmopan, Belize": 17.25,
        "Merida, Mexico": 20.96,
        "Miami, Florida": 25.76,
        "Charlotte, North Carolina": 35.22,
        "State college, Pennsylvania": 40.79,
        "Ottawa, Canada": 45.42
    }
    latitudes = np.array(list(sites.values()))
    
    print("--- Initial State: A Latitudinal Diversity Gradient ---")
    
    # 2. Simulate the initial alpha diversity for each site.
    # We model a scenario where diversity is highest at the equator (~250 species)
    # and decreases with latitude, with some random variation.
    np.random.seed(42) # for reproducible results
    initial_diversity = 250 - 4.5 * latitudes + np.random.normal(0, 10, len(latitudes))
    
    # 3. Calculate the slope of the initial diversity gradient.
    # We use linear regression (polyfit) to find the slope of the line that best 
    # fits the relationship between latitude and diversity.
    initial_slope, initial_intercept = np.polyfit(latitudes, initial_diversity, 1)
    
    print("The initial relationship between diversity and latitude can be modeled.")
    print("The equation is: Diversity = (Slope * Latitude) + Intercept")
    print(f"Initial Equation: Diversity ≈ ({initial_slope:.2f} * Latitude) + {initial_intercept:.2f}")
    print(f"The initial slope is approximately {initial_slope:.2f}, showing a strong negative relationship as expected.\n")

    print("--- Ultimate State: Post-Invasion Community Collapse ---")
    print("After the invasive species out-competes native species everywhere,")
    print("alpha diversity at all sites collapses to a very low, uniform number.")
    
    # 4. Simulate the final alpha diversity.
    # In the ultimate state, the invasive has driven most other species to extinction.
    # We model this by setting the diversity at all sites to a constant low value (e.g., 2).
    final_diversity = np.full(len(sites), 2)
    
    # 5. Calculate the slope of the final diversity gradient.
    final_slope, final_intercept = np.polyfit(latitudes, final_diversity, 1)
    
    print("The final relationship between diversity and latitude is recalculated.")
    print(f"Final Equation: Diversity ≈ ({final_slope:.2f} * Latitude) + {final_intercept:.2f}")
    print(f"The final slope is effectively {final_slope:.2f}, as diversity is now constant across all latitudes.")
    print("\nConclusion: The slope of the diversity gradient will ultimately approach zero.")

solve_diversity_slope()
<<<0>>>