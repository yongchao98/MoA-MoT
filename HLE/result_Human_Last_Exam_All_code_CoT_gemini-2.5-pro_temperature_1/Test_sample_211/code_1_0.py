import numpy as np
from scipy import stats

def solve_diversity_slope():
    """
    Models the change in the latitudinal diversity gradient of insects
    following an invasion, according to Hubbell's Unified Neutral Theory.
    """
    # Plan Step 1: Define the sites and their approximate latitudes (in degrees).
    sites = {
        "Tena, Ecuador": 0,
        "Cali, Colombia": 3,
        "Barro Colorado Island, Panama": 9,
        "Sarapiqui, Costa Rica": 10,
        "Managua, Nicaragua": 12,
        "Yoro, Honduras": 15,
        "Belmopan, Belize": 17,
        "Merida, Mexico": 21,
        "Miami, Florida": 26,
        "Charlotte, North Carolina": 35,
        "State College, Pennsylvania": 41,
        "Ottawa, Canada": 45
    }

    latitudes = np.array(list(sites.values()))
    site_names = list(sites.keys())

    print("--- Initial State: Pre-Invasion ---")
    print("Initially, there is a strong latitudinal diversity gradient (LDG), where species diversity is highest at the equator and decreases towards the poles.\n")

    # Plan Step 2: Model the initial alpha diversity.
    # We'll create a linear relationship: Diversity = Intercept + Slope * Latitude
    # A strong LDG has a steep, negative slope.
    np.random.seed(42) # For reproducible results
    initial_slope_param = -2.8
    initial_intercept_param = 160
    initial_noise = np.random.normal(0, 8, len(latitudes))
    initial_diversity = initial_intercept_param + initial_slope_param * latitudes + initial_noise
    initial_diversity[initial_diversity < 10] = 10 # Ensure no negative/zero diversity

    # Calculate the actual slope from the generated data
    slope_initial, intercept_initial, _, _, _ = stats.linregress(latitudes, initial_diversity)

    print(f"The initial relationship between diversity and latitude can be described by the equation:")
    print(f"Alpha Diversity = ({slope_initial:.2f} * Latitude) + {intercept_initial:.2f}\n")
    print("Initial Alpha Diversity (Number of Species):")
    for i in range(len(site_names)):
        print(f"- {site_names[i]:<30} (Lat: {latitudes[i]:>2}°): {int(initial_diversity[i])}")

    print("\n" + "="*50 + "\n")

    print("--- Ultimate State: Post-Invasion (Hubbell's Neutral Theory) ---")
    print("Under neutral theory, the successful invasive species stochastically replaces native species, leading to biotic homogenization. This reduces diversity everywhere, especially in the species-rich tropics, thus flattening the gradient.\n")

    # Plan Step 3: Model the final alpha diversity after homogenization.
    # The slope becomes much less steep (closer to zero).
    final_slope_param = -0.4
    final_intercept_param = 40
    final_noise = np.random.normal(0, 4, len(latitudes))
    final_diversity = final_intercept_param + final_slope_param * latitudes + final_noise
    final_diversity[final_diversity < 5] = 5 # Ensure a minimum diversity remains

    # Calculate the new slope
    slope_final, intercept_final, _, _, _ = stats.linregress(latitudes, final_diversity)

    print(f"The final relationship between diversity and latitude can be described by the equation:")
    print(f"Alpha Diversity = ({slope_final:.2f} * Latitude) + {intercept_final:.2f}\n")
    print("Final Alpha Diversity (Number of Species):")
    for i in range(len(site_names)):
        print(f"- {site_names[i]:<30} (Lat: {latitudes[i]:>2}°): {int(final_diversity[i])}")

    print("\n" + "="*50 + "\n")

    # Plan Step 4: Conclude by comparing the slopes.
    print("--- Conclusion ---")
    print("By comparing the two scenarios, we can see the effect on the slope of the diversity gradient.")
    print(f"Initial Slope: {slope_initial:.2f}")
    print(f"Final Slope  : {slope_final:.2f}")
    print("\nThe slope of insect diversity across these sites will ultimately decrease, becoming flatter and approaching zero.")

solve_diversity_slope()