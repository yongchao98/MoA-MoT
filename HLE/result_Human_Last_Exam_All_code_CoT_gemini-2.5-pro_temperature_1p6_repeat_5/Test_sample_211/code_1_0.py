import numpy as np
from scipy import stats

def solve():
    """
    Demonstrates the effect of an invasive species on the latitudinal diversity gradient
    under the assumptions of Hubbell's neutral theory.
    """
    # Site names and their approximate latitudes
    sites = {
        "Tena, Ecuador": -0.99,
        "Cali, Colombia": 3.45,
        "Barro Colorado Island, Panama": 9.15,
        "Sarapiqui, Costa Rica": 10.47,
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

    # 1. Create a hypothetical initial diversity gradient (high diversity at low latitudes)
    # We'll model this as: Diversity = 250 - 4 * Latitude
    initial_intercept = 250
    initial_true_slope = -4.0
    # Add some random noise to make it more realistic
    np.random.seed(0)
    initial_diversity = initial_intercept + initial_true_slope * latitudes + np.random.normal(0, 10, len(latitudes))
    initial_diversity[initial_diversity < 10] = 10 # Ensure no negative diversity

    # Calculate the observed initial slope using linear regression
    initial_slope, initial_intercept_fit, _, _, _ = stats.linregress(latitudes, initial_diversity)

    print("--- Initial State (Before Invasion) ---")
    print(f"The initial relationship between diversity and latitude is described by the equation:")
    print(f"Diversity = {initial_slope:.2f} * Latitude + {initial_intercept_fit:.2f}")
    print("\n")


    # 2. Simulate the invasion: Biotic Homogenization
    # The invasive species causes extinctions, reducing diversity more (in absolute terms)
    # where it was highest. We'll model this by reducing diversity by 40% everywhere.
    # This reduces the absolute difference between high and low diversity sites.
    final_diversity = initial_diversity * 0.6

    # Calculate the new slope after homogenization
    final_slope, final_intercept_fit, _, _, _ = stats.linregress(latitudes, final_diversity)

    print("--- Final State (After Invasion) ---")
    print("The invasive species causes extinctions, leading to biotic homogenization.")
    print("The final relationship between diversity and latitude is described by the equation:")
    print(f"Diversity = {final_slope:.2f} * Latitude + {final_intercept_fit:.2f}")
    print("\n")
    
    print("--- Conclusion ---")
    print(f"The slope changed from {initial_slope:.2f} to {final_slope:.2f}.")
    print("The slope has become less negative (closer to zero), indicating that the diversity gradient has become shallower or flatter.")


solve()
<<<The slope of insect diversity will decrease (become less negative/flatter).>>>