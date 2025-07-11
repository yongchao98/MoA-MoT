import numpy as np

def solve_polymer_distance():
    """
    Calculates the predicted distance between polymer chains for ligand 8
    by performing a linear regression on the given data for ligands 1, 2, and 3.
    """
    # Step 1: Define the known data points
    # n = number of carbon atoms in the alkyl chain
    # D = distance between polymer chains in Å
    n_values = np.array([1, 2, 3])
    D_values = np.array([12.4, 12.6, 13.2])

    print("Given data points (Number of Carbons, Distance in Å):")
    for n, D in zip(n_values, D_values):
        print(f"({n}, {D})")
    print("-" * 30)

    # Step 2: Perform linear regression to find the best-fit line D = m*n + c
    # The polyfit function with degree 1 performs a linear regression.
    # It returns the slope (m) and the intercept (c).
    m, c = np.polyfit(n_values, D_values, 1)

    print("Performing linear regression to find the relationship D = m*n + c.")
    print(f"Calculated Slope (m): {m:.4f}")
    print(f"Calculated Intercept (c): {c:.4f}")
    print(f"The best-fit linear equation is: D = {m:.4f} * n + {c:.4f}")
    print("-" * 30)

    # Step 3: Predict the distance for ligand 8 (n=8)
    n_target = 8
    predicted_D = m * n_target + c

    print(f"Predicting the distance for ligand 8 (n = {n_target}).")
    # Print the equation with all numbers as requested
    print(f"Predicted Distance = {m:.4f} * {n_target} + {c:.4f}")
    print(f"Predicted Distance = {m * n_target:.4f} + {c:.4f}")
    print(f"Final Predicted Distance: {predicted_D:.4f} Å")
    print("-" * 30)

    # Step 4: Compare with the answer choices
    answer_choices = {
        'A': 12.5, 'B': 14.0, 'C': 15.7, 'D': 17.0, 'E': 18.5, 'F': 20.0
    }
    
    # Find the closest answer choice
    closest_choice = min(answer_choices, key=lambda k: abs(answer_choices[k] - predicted_D))
    
    print("Comparing the result with the given choices:")
    print(f"The calculated value {predicted_D:.2f} Å is closest to option {closest_choice}: {answer_choices[closest_choice]} Å.")

solve_polymer_distance()