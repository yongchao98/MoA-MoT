import numpy as np

def solve_diversity_slope():
    """
    This function simulates the change in the slope of a latitudinal diversity
    gradient after an invasive species causes biotic homogenization.
    """
    # Define the 12 sites and their approximate latitudes
    sites = {
        "Tena, Ecuador": -1,
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

    # We use the absolute value of latitude as the independent variable,
    # representing distance from the equator.
    latitudes = np.array([abs(lat) for lat in sites.values()])

    # --- Step 1: Model the Initial State ---
    # We model a plausible initial diversity gradient where diversity is highest
    # at the equator (latitude 0) and decreases as latitude increases.
    # Equation form: Diversity = Intercept - Slope * Latitude
    initial_intercept = 150  # Hypothetical max diversity at the equator
    initial_gradient = 2.5   # Diversity decreases by 2.5 species per degree of latitude
    # Add some random noise to make the data more realistic
    np.random.seed(0)
    initial_diversities = initial_intercept - initial_gradient * latitudes + np.random.normal(0, 5, len(latitudes))
    initial_diversities = np.round(initial_diversities).astype(int)

    # Calculate the slope of the best-fit line for the initial data
    m_initial, c_initial = np.polyfit(latitudes, initial_diversities, 1)

    print("--- Initial State: Latitudinal Diversity Gradient ---")
    print(f"The initial calculated slope is: {m_initial:.4f}")
    print("Initial diversity equation representing the gradient:")
    print(f"Diversity = {m_initial:.2f} * |Latitude| + {c_initial:.2f}\n")


    # --- Step 2: Model the Final State after Invasion ---
    # The invasive species causes biotic homogenization. Diversity at all sites
    # drops to a similarly low value, erasing the gradient.
    # We'll model this as a random integer between 2 and 4.
    final_diversities = np.random.randint(2, 5, size=len(latitudes))

    # Calculate the slope of the best-fit line for the final data
    m_final, c_final = np.polyfit(latitudes, final_diversities, 1)

    print("--- Final State: Post-Invasion Homogenization ---")
    print("The final calculated slope approaches zero.")
    print(f"The final calculated slope is: {m_final:.4f}")
    print("Final diversity equation representing the flattened gradient:")
    # We output each number in the final equation as requested.
    print(f"Diversity = {m_final:.2f} * |Latitude| + {c_final:.2f}\n")

    print("Conclusion: The successful invasive species erases the latitudinal diversity gradient,")
    print("causing the slope of diversity across the sites to approach zero.")

solve_diversity_slope()