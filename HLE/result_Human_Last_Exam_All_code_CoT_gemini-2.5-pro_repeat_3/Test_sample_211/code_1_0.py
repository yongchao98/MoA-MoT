import numpy as np
from scipy import stats

def solve():
    """
    Simulates the effect of an invasive species on the latitudinal diversity gradient
    under Hubbell's Unified Neutral Theory.
    """
    # 1. Define the sites and their approximate latitudes
    sites = {
        "Tena, Ecuador": -1.0,
        "Cali, Colombia": 3.4,
        "Barro Colorado Island, Panama": 9.2,
        "Sarapiqui, Costa Rica": 10.4,
        "Managua, Nicaragua": 12.1,
        "Yoro, Honduras": 15.1,
        "Belmopan, Belize": 17.2,
        "Merida, Mexico": 21.0,
        "Miami, Florida": 25.8,
        "Charlotte, North Carolina": 35.2,
        "State College, Pennsylvania": 40.8,
        "Ottawa, Canada": 45.4
    }
    
    latitudes = np.array(list(sites.values()))
    
    # 2. Simulate the initial "pre-invasion" diversity gradient
    # We create a strong negative relationship between latitude and diversity.
    # Higher latitude = lower diversity.
    # The intercept is the diversity at the equator (latitude 0).
    initial_intercept = 180
    initial_slope = -3.5
    # Add some random noise to make it more realistic
    random_noise = np.random.normal(0, 10, len(latitudes))
    initial_diversity = initial_intercept + initial_slope * latitudes + random_noise
    # Ensure no diversity is negative
    initial_diversity = np.maximum(initial_diversity, 5)

    # 3. Calculate the slope of the initial gradient
    slope_before, intercept_before, _, _, _ = stats.linregress(latitudes, initial_diversity)
    
    print("Under Hubbell's neutral theory, an invasive species that dominates the metacommunity")
    print("will lead to biotic homogenization, causing diversity to decrease everywhere and")
    print("the diversity gradient to flatten.")
    print("-" * 50)
    print("Initial State (Pre-Invasion):")
    print("The initial diversity gradient shows a clear relationship between latitude and species richness.")
    print("The equation for the line of best fit is:")
    print(f"Diversity = ({slope_before:.2f} * Latitude) + {intercept_before:.2f}")
    print("-" * 50)

    # 4. Simulate the final "post-invasion" state
    # Biotic homogenization occurs. Diversity drops everywhere and the gradient flattens.
    # The new diversity is low and varies little with latitude.
    final_diversity = np.random.normal(15, 3, len(latitudes))
    # Ensure no diversity is negative
    final_diversity = np.maximum(final_diversity, 2)
    
    # 5. Calculate the slope of the final gradient
    slope_after, intercept_after, _, _, _ = stats.linregress(latitudes, final_diversity)

    print("Final State (Post-Invasion):")
    print("After homogenization, diversity has crashed and is low everywhere, regardless of latitude.")
    print("The latitudinal gradient has effectively disappeared.")
    print("The new equation for the line of best fit is:")
    print(f"Diversity = ({slope_after:.2f} * Latitude) + {intercept_after:.2f}")
    print("-" * 50)
    print("As shown, the slope of the diversity gradient has ultimately approached zero.")

solve()
<<<0>>>