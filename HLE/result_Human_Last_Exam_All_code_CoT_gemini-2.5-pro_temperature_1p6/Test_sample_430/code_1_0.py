import math

def calculate_average_layers():
    """
    Calculates the average number of adsorbed layers per site based on a simplified
    grand canonical ensemble model.

    The model makes the following assumptions based on the problem statement:
    1. Lateral interaction energy (ε_l) is negligible and assumed to be 0.
    2. Adsorption energy for layers n > 1 (ε_n) is assumed to be 0.
    3. The maximum number of layers (k) is not given and is assumed to be 10.
    """

    # Parameters from the problem, in units of k_B*T.
    # βε₁ where β = 1/(k_B*T)
    beta_epsilon_1 = 0.1
    # βμ
    beta_mu = 0.15
    # Assumed maximum number of layers
    k_max = 10

    print(f"Calculation based on the following parameters and assumptions:")
    print(f"β*ε_1 = {beta_epsilon_1}")
    print(f"β*μ = {beta_mu}")
    print(f"Assumed maximum number of layers k = {k_max}\n")

    # The term for a site with n layers (n>0) in the partition function is exp(βε₁ + nβμ)
    # Denominator (ξ): This is the single-site grand partition function.
    # ξ = 1 (for n=0) + Sum_{n=1 to k} [exp(βε₁ + nβμ)]
    partition_function_xi = 1.0
    for n in range(1, k_max + 1):
        exponent = beta_epsilon_1 + n * beta_mu
        partition_function_xi += math.exp(exponent)

    # Numerator for the average layer calculation
    # Num = Sum_{n=1 to k} [n * exp(βε₁ + nβμ)]
    numerator = 0.0
    for n in range(1, k_max + 1):
        exponent = beta_epsilon_1 + n * beta_mu
        numerator += n * math.exp(exponent)

    # Calculate the average number of layers ⟨n⟩
    if partition_function_xi == 0:
        average_layers = 0
    else:
        average_layers = numerator / partition_function_xi

    # Outputting each number in the final equation: ⟨n⟩ = Numerator / ξ
    print("--- Final Equation Components ---")
    print(f"Sum term in numerator: {numerator:.4f}")
    print(f"Grand partition function per site (ξ): {partition_function_xi:.4f}")
    print(f"Final Equation: <n> = {numerator:.4f} / {partition_function_xi:.4f}")
    print("\n--- Result ---")
    print(f"The average number of adsorbed layers per site is: {average_layers:.4f}")
    
    # Return the final numerical answer in the required format
    print(f"\n<<<{average_layers:.4f}>>>")

calculate_average_layers()