import numpy as np
from scipy.stats import linregress

def analyze_diversity_slope():
    """
    Simulates the effect of an invasive species on the latitudinal diversity gradient
    based on Hubbell's unified neutral theory.
    """
    # Site data: Name and approximate latitude
    sites = {
        "Tena, Ecuador": -0.99,
        "Cali, Colombia": 3.45,
        "Barro Colorado Island, Panama": 9.15,
        "Sarapiqui, Costa Rica": 10.40,
        "Managua, Nicaragua": 12.13,
        "Yoro, Honduras": 15.13,
        "Belmopan, Belize": 17.25,
        "Merida, Mexico": 20.96,
        "Miami, Florida": 25.76,
        "Charlotte, North Carolina": 35.22,
        "State College, Pennsylvania": 40.79,
        "Ottawa, Canada": 45.42
    }

    latitudes = np.array(list(sites.values()))
    site_names = list(sites.keys())

    # --- Step 1: Simulate the "Before Invasion" diversity gradient ---
    # We create a hypothetical strong negative correlation between latitude and diversity.
    # Diversity = intercept + slope * latitude + random noise
    np.random.seed(42) # for reproducible results
    initial_slope = -4.5
    initial_intercept = 220
    # Add some noise to make it look more like real data
    noise = np.random.normal(0, 15, size=len(latitudes))
    diversity_before = initial_intercept + initial_slope * latitudes + noise
    # Ensure diversity is a positive integer (as it represents a species count)
    diversity_before = np.maximum(10, diversity_before).astype(int)

    # --- Step 2: Calculate the slope of the initial gradient ---
    slope_before, intercept_before, _, _, _ = linregress(latitudes, diversity_before)

    print("--- Analysis of Insect Diversity vs. Latitude ---")
    print("\nInitial State (Before Invasion):")
    print(f"The initial latitudinal gradient is modeled by the equation:")
    # The user requested each number in the final equation to be printed
    print(f"Diversity = {slope_before:.2f} * Latitude + {intercept_before:.2f}")


    # --- Step 3: Simulate the "After Invasion" state ---
    # The invasive species reduces diversity across all sites. We'll simulate this
    # with a 60% reduction, representing the competitive pressure.
    reduction_factor = 0.40 # (1.0 - 0.60)
    diversity_after = diversity_before * reduction_factor
    # Again, ensure diversity is a positive integer
    diversity_after = np.maximum(1, diversity_after).astype(int)

    # --- Step 4: Calculate the slope of the new, disrupted gradient ---
    slope_after, intercept_after, _, _, _ = linregress(latitudes, diversity_after)

    print("\nFuture State (After Invasion):")
    print("The successful invasive species reduces diversity across all sites.")
    print("The new, disrupted gradient is modeled by the equation:")
    # The user requested each number in the final equation to be printed
    print(f"Diversity = {slope_after:.2f} * Latitude + {intercept_after:.2f}")

    # --- Step 5: Compare the slopes and conclude ---
    print("\n--- Conclusion ---")
    print(f"The initial slope was {slope_before:.2f}.")
    print(f"After the invasion, the slope is {slope_after:.2f}.")
    print("\nThe absolute value of the slope has decreased, meaning the line has become flatter.")
    print("Ultimately, the slope of insect diversity observed across these sites will become shallower.")

if __name__ == '__main__':
    analyze_diversity_slope()
