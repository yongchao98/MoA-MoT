import numpy as np
from scipy import stats

def analyze_diversity_slope():
    """
    Analyzes the change in the slope of a latitudinal diversity gradient
    following a simulated biological invasion, based on neutral theory.
    """
    # Step 1: Define the sites and their approximate latitudes (in degrees).
    sites = {
        "Tena, Ecuador": 0.99,
        "Cali, Colombia": 3.45,
        "Barro Colorado Island, Panama": 9.15,
        "Sarapiqui, Costa Rica": 10.42,
        "Managua, Nicaragua": 12.13,
        "Yoro, Honduras": 15.14,
        "Belmopan, Belize": 17.25,
        "Merida, Mexico": 20.97,
        "Miami, Florida": 25.76,
        "Charlotte, North Carolina": 35.23,
        "State College, Pennsylvania": 40.79,
        "Ottawa, Canada": 45.42,
    }
    
    latitudes = np.array(list(sites.values()))
    
    # Step 2: Create a hypothetical initial alpha diversity representing the LDG.
    # We model diversity as being high at the equator (low lat) and decreasing linearly.
    # Initial diversity = 150 species - 2.8 species per degree of latitude.
    initial_slope_true = -2.8
    initial_intercept_true = 150
    # Add some noise to make it more realistic
    noise = np.random.normal(0, 5, len(latitudes))
    initial_diversity = initial_intercept_true + initial_slope_true * latitudes + noise
    initial_diversity[initial_diversity < 5] = 5 # Ensure no negative diversity

    # Step 3: Calculate the slope of the initial diversity gradient.
    slope_initial, intercept_initial, _, _, _ = stats.linregress(latitudes, initial_diversity)
    
    print("--- Initial State: Latitudinal Diversity Gradient ---")
    print(f"The initial relationship between diversity and latitude is modeled by the equation:")
    # The format specifier '.4f' ensures we print the numbers from the equation.
    print(f"Diversity = {slope_initial:.4f} * Latitude + {intercept_initial:.4f}\n")

    # Step 4: Model the final state after the invasion.
    # Under neutral theory, the invasive species drives diversity down at all sites
    # to a uniformly low level. We'll model this as diversity approaching a small number (e.g., 2).
    final_diversity = np.full_like(latitudes, 2.0)
    # Add a tiny bit of noise to avoid a perfectly flat line, which is more realistic
    final_diversity += np.random.normal(0, 0.1, len(latitudes))
    
    # Step 5: Calculate the slope of the final diversity "gradient".
    slope_final, intercept_final, _, _, _ = stats.linregress(latitudes, final_diversity)
    
    print("--- Final State: After Invasion Homogenization ---")
    print("Under neutral theory, the invasive species homogenizes communities, causing diversity to plummet at all latitudes.")
    print("The final relationship between diversity and latitude is now modeled by the equation:")
    # The format specifier '.4f' ensures we print the numbers from the equation.
    print(f"Diversity = {slope_final:.4f} * Latitude + {intercept_final:.4f}\n")
    
    print("--- Conclusion ---")
    print(f"The initial slope of {slope_initial:.4f} was steep and negative, reflecting a strong diversity gradient.")
    print(f"The final slope of {slope_final:.4f} is nearly flat, indicating the gradient has been erased.")
    print("Ultimately, the slope of insect diversity observed across these sites will approach zero.")

# Run the analysis
analyze_diversity_slope()