import math

def solve_adsorption():
    """
    Calculates the average number of adsorbed layers per site for a
    multi-layer adsorption model under simplifying assumptions.
    """
    # Parameters from the user input, in units of kB*T
    # epsilon_1 / (kB*T)
    e1_kT = 0.1
    # mu / (kB*T)
    mu_kT = 0.15

    # Assumptions based on the problem statement
    # Assume k is related to coordination numbers z_l=4, z_inter=4.
    k = 4
    
    # Based on this k, epsilon_l is negligible: (0.02)^4 is very small.
    # We assume lateral interactions are zero.
    el_kT = 0
    
    # We assume epsilon_2 is also negligible.
    e2_kT = 0
    
    # The problem reduces to a finite-layer BET model.
    # beta * (mu + epsilon_1)
    x1 = math.exp(mu_kT + e1_kT)
    
    # beta * (mu + epsilon_2)
    x2 = math.exp(mu_kT + e2_kT)

    # Calculate the partition function z
    # z = 1 + x1 + x1*x2 + x1*x2^2 + ... + x1*x2^(k-1)
    # z = 1 + x1 * (1 - x2^k) / (1 - x2)
    z_val = 1 + x1 * (1 - math.pow(x2, k)) / (1 - x2)
    
    # Calculate the average number of layers <k>
    # <k> = (1/z) * sum_{m=1 to k} m * w_m
    # where w_m = x1 * x2^(m-1)
    # sum_term = x1 * (1 - (k+1)*x2^k + k*x2^(k+1)) / (1-x2)^2
    numerator_sum = 1 - (k + 1) * math.pow(x2, k) + k * math.pow(x2, k + 1)
    denominator_sum = math.pow(1 - x2, 2)
    sum_term = x1 * (numerator_sum / denominator_sum)
    
    avg_k = sum_term / z_val

    # Print the derivation with numbers
    print("Under the assumption that k=4 and lateral interactions are negligible, the model simplifies to the BET model for finite layers.")
    print("The key parameters are calculated as:")
    print(f"x\u2081 = exp(\u03BC/k\u208BT + \u03B5\u2081/k\u208BT) = exp({mu_kT:.2f} + {e1_kT:.2f}) = {x1:.4f}")
    print(f"x\u2082 = exp(\u03BC/k\u208BT + \u03B5\u2082/k\u208BT) = exp({mu_kT:.2f} + {e2_kT:.2f}) = {x2:.4f}")
    print("\nThe single-site partition function z is:")
    print(f"z = 1 + x\u2081 * (1 - x\u2082\u2074) / (1 - x\u2082) = 1 + {x1:.4f} * (1 - {x2:.4f}\u2074) / (1 - {x2:.4f}) = {z_val:.4f}")
    
    print("\nThe average number of adsorbed layers per site \u27E8k\u27E9 is given by:")
    print("\u27E8k\u27E9 = (1/z) * \u2211\u2098\u208C\u2081\u2074 m\u22C5w\u2098")
    print(f"\u27E8k\u27E9 = (1/z) * x\u2081 * [1 - (k+1)x\u2082\u2075 + kx\u2082\u2075\u207A\u00B9] / (1-x\u2082)\u00B2")
    
    k_plus_1_xk = (k + 1) * math.pow(x2, k)
    k_xk_plus_1 = k * math.pow(x2, k + 1)
    
    print(f"\u27E8k\u27E9 = (1/{z_val:.4f}) * {x1:.4f} * [1 - {k+1:.0f}*({x2:.4f})\u2074 + {k:.0f}*({x2:.4f})\u2075] / (1 - {x2:.4f})\u00B2")
    print(f"\u27E8k\u27E9 = ({x1/z_val:.4f}) * [{1 - k_plus_1_xk:.4f} + {k_xk_plus_1:.4f}] / {math.pow(1 - x2, 2):.4f}")
    print(f"\u27E8k\u27E9 = ({x1/z_val:.4f}) * [{numerator_sum:.4f}] / {denominator_sum:.4f}")
    print(f"\u27E8k\u27E9 = ({x1/z_val:.4f}) * {numerator_sum/denominator_sum:.4f}")
    print(f"\u27E8k\u27E9 = {avg_k:.4f}")
    
    # Final answer in the required format
    print(f"\n<<<{avg_k:.4f}>>>")

solve_adsorption()