import math

def calculate_average_layers():
    """
    Calculates the average number of adsorbed layers per site based on a
    simplified finite-layer BET model.
    """
    # Parameters from the problem statement (in units of k_B*T)
    epsilon_1_tilde = 0.1
    mu_tilde = 0.15
    z_l = 4

    # Assumptions to make the problem solvable
    # 1. Adsorption energy for subsequent layers equals the chemical potential.
    epsilon_s_tilde = mu_tilde
    # 2. Maximum number of layers k is assumed to be 4.
    k_max = 4

    # The lateral interaction energy epsilon_l depends on k_max.
    # epsilon_l_tilde = (0.02)**k_max -> this is very small, so lateral
    # interactions are negligible (w = z_l * epsilon_l_tilde approx 0).
    # We proceed by setting the lateral interaction term to zero.

    # Calculate equilibrium constants
    K1 = math.exp(epsilon_1_tilde + mu_tilde)
    Ks = math.exp(epsilon_s_tilde + mu_tilde)

    # The final equation for the average number of layers <k> is:
    # <k> = sum_{i=1}^{k_max} i * y_i
    # where y_i is the fraction of sites with exactly i layers.
    # y_i can be found from K1 and Ks.

    # Calculate the geometric sum part of the normalization factor
    # sum_{j=0}^{k_max-1} (Ks)^j
    if Ks == 1:
        geom_sum = k_max
    else:
        geom_sum = (1 - Ks**k_max) / (1 - Ks)

    # Calculate y1, the fraction of sites with exactly one layer
    y1 = K1 / (1 + K1 * geom_sum)

    # Calculate y_i for all i from 1 to k_max
    y = [0.0] * k_max
    y[0] = y1
    for i in range(1, k_max):
        y[i] = y[i-1] * Ks
    
    # The final equation is <k> = y_1 * 1 + y_2 * 2 + y_3 * 3 + y_4 * 4
    # Calculate the average number of layers
    avg_k = 0
    equation_str_parts = []
    for i in range(k_max):
        avg_k += (i + 1) * y[i]
        equation_str_parts.append(f"({i+1}) * {y[i]:.4f}")

    print("--- Calculation Details ---")
    print(f"Assumed maximum layers, k = {k_max}")
    print(f"Assumed subsequent layer energy, epsilon_s = {epsilon_s_tilde:.2f} k_B*T")
    print(f"Equilibrium constant for the first layer, K1 = {K1:.4f}")
    print(f"Equilibrium constant for subsequent layers, Ks = {Ks:.4f}")
    print("\nFractions of sites with exactly i layers (y_i):")
    for i in range(k_max):
        print(f"  y_{i+1} = {y[i]:.4f}")
    
    print("\nFinal equation for the average number of layers <k>:")
    final_equation = " + ".join(equation_str_parts)
    print(f"<k> = {final_equation}")
    print(f"<k> = {avg_k}")
    print("--------------------------")
    
    # Return the final numerical answer as requested by the user prompt
    return avg_k

# Run the calculation and print the final result in the required format
final_answer = calculate_average_layers()
# print(f"\n<<< {final_answer} >>>")
# The problem asks for the answer to be at the very end.
# The printing inside the function provides the step-by-step thinking process requested.

if __name__ == '__main__':
    pass # The value is captured in final_answer for the last print statement.