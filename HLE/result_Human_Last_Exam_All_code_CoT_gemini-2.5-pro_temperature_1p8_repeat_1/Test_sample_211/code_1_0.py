import numpy as np

def calculate_and_print_regression(latitudes, diversity_data, state_label):
    """Calculates and prints the slope and equation for a given dataset."""
    
    # Calculate the slope and intercept using a linear fit (degree 1 polynomial)
    slope, intercept = np.polyfit(latitudes, diversity_data, 1)

    print(f"--- {state_label} ---")
    if state_label == "Initial State":
        print("The initial relationship between diversity and latitude shows a strong negative slope.")
    else:
        print("After homogenization, diversity is uniformly low, and the latitudinal gradient is erased.")

    print(f"Calculated Slope (m): {slope:.4f}")
    print(f"Calculated Intercept (b): {intercept:.4f}")
    
    # Print the full equation: Diversity = m * Latitude + b
    # We output each number in the equation as requested.
    print(f"Equation: Diversity = {slope:.4f} * Latitude + {intercept:.4f}")
    print("\n" + "-"*25 + "\n")

# Define the study sites and their approximate latitudes
sites_latitudes = {
    "Tena, Ecuador": -0.99,
    "Cali, Colombia": 3.45,
    "Barro Colorado Island, Panama": 9.15,
    "Sarapiqui, Costa Rica": 10.4,
    "Managua, Nicaragua": 12.13,
    "Yoro, Honduras": 15.1,
    "Belmopan, Belize": 17.25,
    "Merida, Mexico": 20.96,
    "Miami, Florida": 25.76,
    "Charlotte, North Carolina": 35.22,
    "State College, Pennsylvania": 40.79,
    "Ottawa, Canada": 45.42
}

latitudes_array = np.array(list(sites_latitudes.values()))

# 1. Initial State: A clear Latitudinal Diversity Gradient (LDG)
# Higher diversity (more species) at lower latitudes.
initial_diversity_array = np.array([305, 280, 250, 245, 240, 220, 210, 190, 170, 120, 100, 80])
calculate_and_print_regression(latitudes_array, initial_diversity_array, "Initial State")

# 2. Final State: After the invasive species has dominated all sites under Neutral Theory
# Alpha diversity collapses to ~1 everywhere, erasing the gradient.
# The small variations represent stochastic noise.
final_diversity_array = np.array([1.2, 1.1, 1.3, 1.2, 1.1, 1.2, 1.3, 1.1, 1.2, 1.1, 1.3, 1.2])
calculate_and_print_regression(latitudes_array, final_diversity_array, "Ultimate Final State")

print("As shown by the calculation for the final state, the slope ultimately approaches 0.")
