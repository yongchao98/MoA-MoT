import numpy as np
from scipy import stats

def analyze_diversity_gradient():
    """
    Simulates and analyzes the change in the latitudinal diversity gradient
    due to a widespread invasive species, according to Hubbell's Neutral Theory.
    """
    # Step 1: Define the study sites and their approximate latitudes.
    # Latitudes increase as we move north from the equator.
    sites = [
        "Tena, Ecuador", "Cali, Colombia", "Barro Colorado Island, Panama", 
        "Sarapiqui, Costa Rica", "Managua, Nicaragua", "Yoro, Honduras", 
        "Belmopan, Belize", "Merida, Mexico", "Miami, Florida", 
        "Charlotte, North Carolina", "State college, Pennsylvania", "Ottawa, Canada"
    ]
    latitudes = np.array([
        -1.0, 3.4, 9.1, 10.4, 12.1, 15.1, 
        17.2, 21.0, 25.8, 35.2, 40.8, 45.4
    ])

    # Step 2: Simulate the initial (pre-invasion) diversity gradient.
    # We create a linear model: Diversity = intercept - slope * latitude
    # This shows high diversity at low latitudes (near equator).
    # We add some random noise to make it more realistic.
    np.random.seed(0) # for reproducibility
    initial_intercept = 250
    initial_slope_param = -4.0
    initial_diversity = initial_intercept + initial_slope_param * latitudes + np.random.normal(0, 15, len(latitudes))
    initial_diversity[initial_diversity < 10] = 10 # ensure no negative diversity

    # Step 3: Analyze the initial gradient by calculating its slope.
    slope_before, intercept_before, r_val_before, _, _ = stats.linregress(latitudes, initial_diversity)
    
    print("--- Initial State (Pre-Invasion) ---")
    print("A simulation of the initial latitudinal diversity gradient.")
    print("Diversity is highest near the equator (low latitude) and decreases northward.")
    print("The equation for the trendline is: Diversity = (Slope * Latitude) + Intercept")
    print(f"Initial Equation: Diversity = ({slope_before:.4f} * Latitude) + {intercept_before:.4f}")
    print(f"The initial R-squared value is: {r_val_before**2:.4f}\n")

    # Step 4: Simulate the ultimate outcome under Neutral Theory (post-invasion).
    # The invasive species causes homogenization, making diversity values across
    # sites converge toward their average. We simulate this by shifting each
    # diversity value closer to the mean.
    mean_diversity = np.mean(initial_diversity)
    homogenization_factor = 0.85 # Strong convergence (85%)
    final_diversity = initial_diversity * (1 - homogenization_factor) + mean_diversity * homogenization_factor

    # Step 5: Analyze the final gradient.
    slope_after, intercept_after, r_val_after, _, _ = stats.linregress(latitudes, final_diversity)
    
    print("--- Final State (Post-Invasion) ---")
    print("Under Neutral Theory, the invasive species homogenizes communities over time.")
    print("Diversity levels across all sites become more similar, weakening the gradient.")
    print(f"Final Equation: Diversity = ({slope_after:.4f} * Latitude) + {intercept_after:.4f}")
    print(f"The final R-squared value is: {r_val_after**2:.4f}\n")
    
    print("--- Conclusion ---")
    print(f"The absolute value of the slope has decreased from {abs(slope_before):.4f} to {abs(slope_after):.4f}.")
    print("This indicates the slope of insect diversity has flattened, moving closer to zero.")

analyze_diversity_gradient()