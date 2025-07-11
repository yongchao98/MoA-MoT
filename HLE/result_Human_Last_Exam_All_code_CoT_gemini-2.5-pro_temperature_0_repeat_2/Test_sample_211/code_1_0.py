import numpy as np
from scipy import stats

def analyze_diversity_slope():
    """
    Analyzes the change in the slope of a latitudinal diversity gradient
    following a simulated invasive species introduction, based on the principles
    of Hubbell's unified neutral theory.
    """
    # Site names and their approximate latitudes
    sites = {
        "Tena, Ecuador": -0.99,
        "Cali, Colombia": 3.44,
        "Barro Colorado Island, Panama": 9.15,
        "Sarapiqui, Costa Rica": 10.43,
        "Managua, Nicaragua": 12.13,
        "Yoro, Honduras": 15.09,
        "Belmopan, Belize": 17.25,
        "Merida, Mexico": 20.97,
        "Miami, Florida": 25.76,
        "Charlotte, North Carolina": 35.22,
        "State college, Pennsylvania": 40.79,
        "Ottawa, Canada": 45.42
    }

    latitudes = np.array(list(sites.values()))

    # 1. Initial State: A strong latitudinal diversity gradient.
    # Let's assume a hypothetical diversity that is high at the equator and decreases with latitude.
    # Equation: Diversity = 250 - 4 * |Latitude|
    initial_diversity = 250 - 4 * np.abs(latitudes)

    # Calculate the initial slope using linear regression
    slope_initial, intercept_initial, _, _, _ = stats.linregress(np.abs(latitudes), initial_diversity)

    print("--- Initial State Before Invasion ---")
    print("The initial diversity gradient shows higher diversity near the equator.")
    print("The equation for the diversity gradient is approximately:")
    print(f"Diversity = {intercept_initial:.2f} + ({slope_initial:.2f} * |Latitude|)\n")


    # 2. Post-Invasion State: Biotic Homogenization
    # The invasive species reduces diversity at all sites, causing the gradient to flatten.
    # The range of diversity values shrinks.
    # New Equation: Diversity = 180 - 2.5 * |Latitude|
    post_invasion_diversity = 180 - 2.5 * np.abs(latitudes)

    # Calculate the new slope
    slope_final, intercept_final, _, _, _ = stats.linregress(np.abs(latitudes), post_invasion_diversity)

    print("--- Post-Invasion State (Hubbell's Neutral Theory) ---")
    print("After the invasive species homogenizes communities, the gradient flattens.")
    print("The equation for the new diversity gradient is approximately:")
    print(f"Diversity = {intercept_final:.2f} + ({slope_final:.2f} * |Latitude|)\n")

    # 3. Conclusion
    print("--- Conclusion ---")
    print(f"The initial slope was {slope_initial:.2f}.")
    print(f"The final slope is {slope_final:.2f}.")
    print("The slope has become less negative (it increased from -4.00 to -2.50), moving closer to zero.")
    print("This indicates that the diversity gradient has become flatter or less steep.")

# Run the analysis
analyze_diversity_slope()