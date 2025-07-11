import numpy as np

def calculate_and_print_model(sites, latitudes, diversities, model_name):
    """
    Calculates the linear regression for the given data and prints the results.
    """
    # Perform linear regression (degree 1 polynomial fit)
    # It returns [slope, intercept]
    slope, intercept = np.polyfit(latitudes, diversities, 1)

    print(f"--- {model_name} State ---")
    print("The relationship between insect diversity and latitude can be modeled by the equation:")
    # We print the numbers for the equation as requested
    print(f"Diversity = {slope:.2f} * |Latitude| + {intercept:.2f}\n")
    return slope

def main():
    """
    Main function to run the simulation.
    """
    # Step 1: Establish study sites and their approximate absolute latitudes
    sites = {
        "Tena, Ecuador": 0.99,
        "Cali, Colombia": 3.45,
        "Barro Colorado Island, Panama": 9.15,
        "Sarapiqui, Costa Rica": 10.40,
        "Managua, Nicaragua": 12.13,
        "Yoro, Honduras": 15.14,
        "Belmopan, Belize": 17.25,
        "Merida, Mexico": 20.97,
        "Miami, Florida": 25.76,
        "Charlotte, North Carolina": 35.23,
        "State college, Pennsylvania": 40.79,
        "Ottawa, Canada": 45.42
    }
    site_names = list(sites.keys())
    # Use absolute latitude as we are measuring distance from the equator
    abs_latitudes = np.array(list(sites.values()))

    # Step 2: Simulate the Initial "Pre-Invasion" State
    # We model a strong latitudinal diversity gradient.
    # High diversity at the equator (intercept) that drops off quickly (steep slope).
    # Adding a little noise to make the regression more realistic.
    np.random.seed(0)
    initial_intercept = 250
    initial_slope_param = -4.5
    noise = np.random.normal(0, 10, len(abs_latitudes))
    initial_diversities = initial_intercept + initial_slope_param * abs_latitudes + noise

    # Step 3: Simulate the Final "Post-Invasion" State
    # Biotic homogenization reduces diversity more in the tropics, flattening the gradient.
    # The intercept is lower and the slope is less steep (closer to zero).
    final_intercept = 120
    final_slope_param = -1.5
    noise = np.random.normal(0, 5, len(abs_latitudes))
    final_diversities = final_intercept + final_slope_param * abs_latitudes + noise

    # Step 4: Calculate and Compare Slopes
    print("Modeling the change in the latitudinal diversity gradient after a widespread invasion.\n")

    initial_slope = calculate_and_print_model(site_names, abs_latitudes, initial_diversities, "Initial (Pre-Invasion)")
    final_slope = calculate_and_print_model(site_names, abs_latitudes, final_diversities, "Final (Post-Invasion)")

    # Step 5: Conclusion
    print("--- Conclusion ---")
    print(f"The initial slope was {initial_slope:.2f}, showing a steep diversity gradient.")
    print(f"After the invasion, the final slope is {final_slope:.2f}.")
    print("\nThe invasion causes biotic homogenization, disproportionately reducing diversity in the species-rich tropics.")
    print("This flattens the latitudinal diversity gradient, causing the slope to become less steep (i.e., less negative, increasing towards zero).")


if __name__ == "__main__":
    main()

<<<The slope of the latitudinal diversity gradient will become less steep, increasing in value as it approaches zero.>>>