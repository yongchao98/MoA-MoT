import numpy as np

def calculate_average_layers():
    """
    Calculates the average number of adsorbed layers per site based on a
    simplified model derived from the problem description.

    The model makes the following assumptions due to ambiguities in the prompt:
    1.  The maximum number of layers, k_max, is assumed to be equal to the
        lateral coordination number, z_l. So, k_max = 4.
    2.  The lateral interaction energy, ε_l = (0.02)^k_max * k_B*T, is negligibly
        small. Thus, we assume ε_l = 0 and no mean-field effects.
    3.  The vertical interaction energy for layers j>1, ε_v, is also assumed to
        be negligible (ε_v = 0).
    4.  The energy of a site with j particles is E_j = -ε_1 for j>=1, and E_0 = 0.
    5.  All parameters given in the prompt are used, including the positive
        chemical potential μ = 0.15 * k_B*T.
    """
    # Parameters from the problem or our assumptions
    k_B_T_Joules = 1.380649e-23 * 318
    eps1_over_kT = 0.1
    mu_over_kT = 0.15
    k_max = 4 # Assumption: k_max = z_l

    print("--- Model Parameters ---")
    print(f"Maximum number of layers (k_max): {k_max}")
    print(f"Monolayer interaction energy (ε₁/k_B·T): {eps1_over_kT}")
    print(f"Chemical potential (μ/k_B·T): {mu_over_kT}")
    print(f"Higher layer interaction energy (ε_v/k_B·T): 0.0 (assumed negligible)")
    print("-" * 26)

    # Statistical weights w_j = exp(-(E_j - j*μ)/(k_B*T))
    # E_j = -ε_1 for j>=1, E_0 = 0
    # So, w_j = exp((ε_1 + j*μ)/(k_B*T)) for j>=1, and w_0 = 1
    
    weights = np.zeros(k_max + 1)
    weights[0] = 1.0
    for j in range(1, k_max + 1):
        exponent = eps1_over_kT + j * mu_over_kT
        weights[j] = np.exp(exponent)

    # Calculate the partition function (denominator)
    z_site = np.sum(weights)

    # Calculate the weighted sum of layers (numerator)
    numerator = np.sum(j * weights[j] for j in range(k_max + 1))
    
    # Calculate the average number of layers
    average_k = numerator / z_site
    
    # --- Output the final equation and result ---
    print("The average number of layers per site, <k>, is calculated by:")
    print("<k> = (Σ_{j=0}^{k_max} j * w_j) / (Σ_{j=0}^{k_max} w_j)")
    print("where the statistical weights w_j are:")
    print("w_0 = 1")
    for j in range(1, k_max + 1):
        print(f"w_{j} = exp((ε₁/k_B·T) + {j}*(μ/k_B·T)) = exp(({eps1_over_kT}) + {j}*({mu_over_kT})) = {weights[j]:.4f}")
    
    print("\n--- Calculation ---")
    print(f"Numerator = (0*w_0 + 1*w_1 + 2*w_2 + 3*w_3 + 4*w_4)")
    print(f"Numerator = (0*{weights[0]:.4f} + 1*{weights[1]:.4f} + 2*{weights[2]:.4f} + 3*{weights[3]:.4f} + 4*{weights[4]:.4f}) = {numerator:.4f}")

    print(f"Denominator = (w_0 + w_1 + w_2 + w_3 + w_4)")
    print(f"Denominator = ({weights[0]:.4f} + {weights[1]:.4f} + {weights[2]:.4f} + {weights[3]:.4f} + {weights[4]:.4f}) = {z_site:.4f}")

    print("\n--- Final Result ---")
    print(f"The average number of adsorbed layers per site is <k> = {numerator:.4f} / {z_site:.4f} = {average_k:.4f}")
    
    # The final answer format for the platform
    print(f"\n<<<{average_k:.4f}>>>")

if __name__ == '__main__':
    calculate_average_layers()