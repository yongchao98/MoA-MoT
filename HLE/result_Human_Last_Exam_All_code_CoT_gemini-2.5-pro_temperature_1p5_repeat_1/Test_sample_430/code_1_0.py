import numpy as np

def calculate_average_layers():
    """
    Calculates the average number of adsorbed layers per site based on the derived model.
    """
    # Parameters from the problem (in units of k_B*T)
    e1_over_kT = 0.1
    mu_over_kT = 0.15

    # Assumptions based on the plan
    k_max = 4
    # For layers j > 1, epsilon_j = mu, so (epsilon_j - mu) / kT = 0
    ej_over_kT_j_gt_1 = mu_over_kT

    print("--- Model Assumptions and Setup ---")
    print(f"Maximum number of layers (k): {k_max}")
    print("Lateral interactions (epsilon_l) are assumed to be negligible.")
    print(f"epsilon_1 / (k_B*T) = {e1_over_kT}")
    print(f"mu / (k_B*T) = {mu_over_kT}")
    print(f"For layers j > 1, epsilon_j / (k_B*T) is assumed to be equal to mu / (k_B*T) = {ej_over_kT_j_gt_1}")
    print("-" * 35)

    # Calculate the statistical weight terms X_j
    # x_m = exp(-(epsilon_m/kT - mu/kT))
    x_terms = []
    # For layer m=1
    x1 = np.exp(-(e1_over_kT - mu_over_kT))
    x_terms.append(x1)
    # For layers m > 1
    for _ in range(2, k_max + 1):
        x_m = np.exp(-(ej_over_kT_j_gt_1 - mu_over_kT))
        x_terms.append(x_m)
    
    # X_j = product of x_m from m=1 to j
    X_j_terms = [1.0] # Corresponds to j=0
    current_product = 1.0
    for j in range(k_max):
        current_product *= x_terms[j]
        X_j_terms.append(current_product)

    # Calculate the partition function z_site
    z_site = sum(X_j_terms)
    
    # Calculate the numerator for the average number of layers
    numerator = 0
    for j in range(1, k_max + 1):
        numerator += j * X_j_terms[j]

    # Calculate the average number of layers
    avg_layers = numerator / z_site
    
    # Print the equation and the components
    print("--- Calculation of Average Layers <h> ---")
    print("The formula is: <h> = Sum(j * X_j) / Sum(X_j), for j=0 to k")
    print("\nCalculated statistical weights X_j = exp(-beta * (E_j - j*mu)):")
    for j, Xj in enumerate(X_j_terms):
        print(f"  X_{j}: {Xj:.5f}")

    print(f"\nThe equation for the average number of layers is:")
    numerator_str = " + ".join([f"{j}*({X_j_terms[j]:.5f})" for j in range(1, k_max + 1)])
    denominator_str = " + ".join([f"{X_j_terms[j]:.5f}" for j in range(k_max + 1)])
    print(f"<h> = ({numerator_str}) / ({denominator_str})")

    print(f"\nNumerator = {numerator:.5f}")
    print(f"Denominator (z_site) = {z_site:.5f}")
    
    print("\n--- Final Result ---")
    print(f"The average number of adsorbed layers per site is: {avg_layers:.4f}")
    
    return avg_layers

# Run the calculation
final_answer = calculate_average_layers()
print(f"<<<{final_answer:.4f}>>>")
