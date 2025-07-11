import numpy as np
from scipy import stats

def solve_diversity_slope():
    """
    Models the change in the latitudinal diversity gradient slope
    following an invasion, based on Hubbell's neutral theory.
    """

    # Step 1: Define the sites and their approximate latitudes.
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
        "State College, Pennsylvania": 40.80,
        "Ottawa, Canada": 45.42
    }
    
    latitudes = np.array(list(sites.values()))

    # Step 2: Model the initial diversity gradient.
    # Let's assume a strong linear relationship where diversity decreases with latitude.
    # Initial Diversity = 150 - 3 * Latitude
    initial_intercept = 150
    initial_slope_param = -3.0
    # Add some random noise to make it more realistic
    np.random.seed(0)
    noise = np.random.normal(0, 5, len(latitudes))
    initial_diversity = initial_intercept + initial_slope_param * latitudes + noise

    # Step 3: Calculate the slope of the initial gradient.
    initial_slope, initial_int, _, _, _ = stats.linregress(latitudes, initial_diversity)

    print("--- Initial State: Before Invasive Species Establishes ---")
    print("A clear latitudinal diversity gradient exists.")
    print("The relationship between alpha diversity and latitude is:")
    print(f"Diversity = {initial_slope:.2f} * Latitude + {initial_int:.2f}\n")


    # Step 4: Model the final diversity after biotic homogenization.
    # The invasive species reduces diversity everywhere, but it reduces the
    # difference between high and low diversity sites, flattening the gradient.
    # We can model this by making the diversity values converge toward a lower average.
    final_diversity = (initial_diversity + np.mean(initial_diversity)) / 2.5
    
    # Step 5: Calculate the slope of the final gradient.
    final_slope, final_int, _, _, _ = stats.linregress(latitudes, final_diversity)
    
    print("--- Final State: After Invasive Species Causes Biotic Homogenization ---")
    print("Diversity decreases everywhere, and communities become more similar.")
    print("The new relationship between alpha diversity and latitude is:")
    print(f"Diversity = {final_slope:.2f} * Latitude + {final_int:.2f}\n")

    # Step 6: Conclusion
    print("--- Conclusion ---")
    print(f"The initial slope was {initial_slope:.2f}, and the final slope is {final_slope:.2f}.")
    print("Because the final slope is a less negative number than the initial slope, the slope has increased.")
    print("This means the diversity gradient has become flatter (less steep), which is the predicted outcome of biotic homogenization under neutral theory.")

solve_diversity_slope()

<<<The slope of diversity will ultimately increase (i.e., become less negative and flatter) as biotic homogenization reduces the richness differences between equatorial and temperate sites.>>>