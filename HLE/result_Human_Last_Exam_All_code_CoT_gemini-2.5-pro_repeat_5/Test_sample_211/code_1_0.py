import numpy as np

def solve_diversity_gradient():
    """
    Simulates the effect of an invasive species on the latitudinal diversity gradient
    based on Hubbell's unified neutral theory.
    """
    # Define the sites and their approximate absolute latitudes
    sites = [
        "Tena, Ecuador", "Cali, Colombia", "Barro Colorado Island, Panama",
        "Sarapiqui, Costa Rica", "Managua, Nicaragua", "Yoro, Honduras",
        "Belmopan, Belize", "Merida, Mexico", "Miami, Florida",
        "Charlotte, North Carolina", "State college, Pennsylvania", "Ottawa, Canada"
    ]
    # Using absolute latitude, as the gradient is symmetrical around the equator
    abs_latitudes = np.array([
        0.99, 3.45, 9.15, 10.4, 12.13, 15.14, 17.25, 20.97, 25.76, 35.23, 40.80, 45.42
    ])

    # 1. Model the initial diversity gradient.
    # We assume a linear relationship: Diversity = Intercept + Slope * Latitude
    # A strong gradient has a high intercept and a steep negative slope.
    np.random.seed(42) # for reproducible results
    initial_intercept = 160
    initial_slope_parameter = -2.5
    # Add some random noise to make it more realistic
    initial_diversity = initial_intercept + initial_slope_parameter * abs_latitudes + np.random.normal(0, 5, len(abs_latitudes))
    initial_diversity = np.round(initial_diversity).astype(int)

    # 2. Model the invasion causing biotic homogenization.
    # The number of species is reduced proportionally across all sites.
    # Let's assume a 30% reduction in local species richness due to the invasive.
    reduction_factor = 0.70
    final_diversity = np.round(initial_diversity * reduction_factor).astype(int)

    # 3. Calculate the actual slope of the gradient before and after the invasion
    # using linear regression (numpy's polyfit with degree=1).
    initial_slope, initial_intercept_fit = np.polyfit(abs_latitudes, initial_diversity, 1)
    final_slope, final_intercept_fit = np.polyfit(abs_latitudes, final_diversity, 1)

    # 4. Print the explanation and results.
    print("--- Analysis based on Hubbell's Unified Neutral Theory ---")
    print("The invasive insect will cause 'biotic homogenization,' reducing diversity at all sites.")
    print("This reduction is proportionally larger in absolute terms at the high-diversity equatorial sites.")
    print("This causes the latitudinal diversity gradient to flatten.\n")
    
    print("--- Simulation Results ---")
    print(f"Initial Diversity (from equator to pole): {list(initial_diversity)}")
    print(f"Final Diversity (after invasion):         {list(final_diversity)}\n")

    print("The equation for the diversity gradient is 'Diversity = Intercept + (Slope * Absolute_Latitude)'.")
    
    print("\nInitial Gradient Equation:")
    # Printing each number in the final equation as requested
    print(f"Diversity = {initial_intercept_fit:.2f} + ({initial_slope:.2f} * Absolute_Latitude)")

    print("\nFinal Gradient Equation (after invasion):")
    # Printing each number in the final equation as requested
    print(f"Diversity = {final_intercept_fit:.2f} + ({final_slope:.2f} * Absolute_Latitude)")

    print(f"\nConclusion: The slope became less steep, changing from {initial_slope:.2f} to {final_slope:.2f}, which is closer to zero.")

solve_diversity_gradient()
<<<The slope of the diversity gradient will become less steep, approaching zero.>>>