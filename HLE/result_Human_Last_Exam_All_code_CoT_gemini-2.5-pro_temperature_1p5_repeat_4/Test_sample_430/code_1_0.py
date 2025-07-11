import math

def solve_adsorption():
    """
    Calculates the average number of adsorbed layers per site
    based on the provided model and parameters.
    """
    # Parameters from the problem statement, in units of kB*T
    epsilon_1_norm = 0.1
    mu_norm = 0.15
    z_l = 4
    z_inter = 4

    # --- Step-by-step reasoning for model simplification ---
    # 1. Assume k is related to the given coordination numbers.
    k = 4
    
    # 2. Check lateral interaction energy.
    # epsilon_l_norm = (0.02)**k which is (0.02)**4 = 1.6e-7.
    # The lateral interaction term z_l * epsilon_l_norm is ~6.4e-7,
    # which is negligible compared to other energy scales (0.1, 0.15).
    # So we neglect lateral interactions.
    
    # 3. Assume the saturation condition for upper layers (a standard assumption
    # in BET-like models), so the energy of subsequent layers equals
    # the chemical potential.
    # epsilon_s = mu, which means x_s = exp(beta*(mu - epsilon_s)) = exp(0) = 1.

    # --- Calculation ---
    # The simplified formula for the average number of layers <n> is:
    # <n> = (x1 * k * (k+1) / 2) / (1 + k * x1)
    # where x1 = exp(beta * (mu - epsilon_1)) = exp(mu_norm - epsilon_1_norm)

    # Calculate x1
    x1 = math.exp(mu_norm - epsilon_1_norm)

    # Calculate the numerator and denominator terms
    sum_j = k * (k + 1) / 2
    numerator = x1 * sum_j
    denominator = 1 + k * x1

    # Calculate the average number of layers
    average_layers = numerator / denominator

    # Output the final equation with all numbers
    print("Based on the assumptions:")
    print(f"  Maximum number of layers, k = {k}")
    print(f"  Energy of second and further layers, epsilon_s = mu")
    print(f"  Lateral interactions are negligible.")
    print("\nThe calculation follows the formula: <k> = (x1 * (k * (k+1) / 2)) / (1 + k * x1)")
    print("\nIntermediate values:")
    print(f"  k = {k}")
    print(f"  x1 = exp({mu_norm} - {epsilon_1_norm}) = {x1:.5f}")
    print(f"  k * (k+1) / 2 = {sum_j}")
    print("\nFinal Equation:")
    print(f"<k> = ({x1:.5f} * {sum_j}) / (1 + {k} * {x1:.5f})")
    print(f"    = {numerator:.5f} / {denominator:.5f}")
    print(f"    = {average_layers:.5f}")
    print("\n<<<" + f"{average_layers}" + ">>>")

solve_adsorption()