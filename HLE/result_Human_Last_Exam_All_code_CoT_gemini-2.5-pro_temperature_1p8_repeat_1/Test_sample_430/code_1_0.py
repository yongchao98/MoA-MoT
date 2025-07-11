import math

def calculate_adsorption():
    """
    Calculates the average number of adsorbed layers per site based on the
    provided physical model and parameters.
    """
    # 1. Define normalized parameters from the problem description.
    # We assume k_max (maximum number of layers) = 2, as it's not specified.
    k_max = 2

    mu_norm = 0.15  # mu / (kB * T)
    e1_norm = 0.1   # epsilon_1 / (kB * T)
    zl = 4          # Lateral coordination number
    z_inter = 4     # Vertical coordination number

    # 2. Calculate normalized interaction energies based on assumptions.
    # Assumption: el_norm depends on k_max
    el_norm = (0.02)**k_max

    # Assumption: e2_norm is derived from z_inter and el_norm
    e2_norm = z_inter * el_norm

    # 3. Calculate effective binding energies including the mean-field term.
    # The mean-field energy correction per particle is zl * el_norm / 2.
    mf_correction = (zl * el_norm) / 2
    e1_eff_norm = e1_norm + mf_correction
    e2_eff_norm = e2_norm + mf_correction

    # 4. Calculate the statistical weight factors x and y.
    # x corresponds to adding the first layer.
    # y corresponds to adding a subsequent layer.
    x = math.exp(mu_norm + e1_eff_norm)
    y = math.exp(mu_norm + e2_eff_norm)

    # 5. Calculate the average number of layers for k_max = 2.
    # The partition function for a single site is q = 1 + x + x*y.
    # The numerator for the average is sum(m * P_m) = 1*x + 2*x*y.
    q = 1 + x + x * y
    numerator = x + 2 * x * y
    average_layers = numerator / q

    # Output the final equation and result
    print(f"The calculation is for a maximum of k={k_max} layers.")
    print("The formula for the average number of layers is: <k> = (x + 2*x*y) / (1 + x + x*y)")
    print(f"Calculated statistical weight for the first layer (x): {x}")
    print(f"Calculated statistical weight for the second layer (y): {y}")
    print("\nFinal Result:")
    print(f"The average number of adsorbed layers per site is: {average_layers}")

    return average_layers

# Execute the calculation and store the result
result = calculate_adsorption()

# Print the final answer in the required format
# Note: The output above already explains the calculation.
# The following line provides the single numerical answer as requested.
# To avoid redundancy, typically you would just have the final line.
# print(f'<<<{result}>>>')
