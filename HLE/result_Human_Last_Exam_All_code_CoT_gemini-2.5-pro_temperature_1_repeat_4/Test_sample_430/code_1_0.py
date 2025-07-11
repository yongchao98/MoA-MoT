import math

def solve_adsorption():
    """
    Calculates the average number of adsorbed layers per site based on the
    grand canonical ensemble and a set of simplifying assumptions derived from
    the problem statement.
    """
    # Parameters given in the problem, in units of k_B*T
    e1_kT = 0.1
    mu_kT = 0.15

    # --- Plan ---
    # 1. We determine the maximum number of layers, k. The problem gives
    #    coordination numbers z_l=4 and z_inter=4. This suggests that the
    #    intended maximum number of layers is k=4.
    k = 4

    # 2. We analyze the interaction energies. The lateral interaction energy
    #    epsilon_l = (0.02)^k * kT is negligible for k=4. We also assume that
    #    for layers i>=2, the adsorption energy equals the chemical potential,
    #    epsilon_i = mu, analogous to adsorption from a saturated phase.

    # 3. We calculate the average number of layers <k> using the formula derived
    #    from the single-site grand partition function under these assumptions:
    #    <k> = [x * k * (k + 1) / 2] / [1 + k * x]
    #    where x = exp(-(e1 - mu)/kT).

    # --- Calculation ---
    print("--- Calculation Steps ---")
    print(f"Based on the problem statement, we assume the maximum number of layers is k = {k}.")

    # The argument of the exponential factor x
    exp_arg = -(e1_kT - mu_kT)
    # The factor x
    x = math.exp(exp_arg)

    # The sum of indices j from 1 to k
    sum_j = k * (k + 1) / 2

    # Numerator of the final expression for <k>
    numerator = x * sum_j

    # Denominator of the final expression for <k>
    denominator = 1 + k * x

    # Final result
    average_layers = numerator / denominator

    print("\nThe final equation for the average number of layers <k> is:")
    print("<k> = (k * (k + 1) / 2) * x / (1 + k * x)")
    print("where x = exp(-(e1/kT - mu/kT))")
    print("\nPlugging in the numbers:")
    print(f"k = {k}")
    print(f"e1/kT = {e1_kT}")
    print(f"mu/kT = {mu_kT}")
    print(f"The exponential term x = exp(-({e1_kT} - {mu_kT})) = exp({exp_arg:.2f}) = {x:.6f}")
    print(f"The sum of indices j from 1 to k = k*(k+1)/2 = {sum_j}")
    print(f"The numerator of the expression for <k> = x * (k*(k+1)/2) = {x:.6f} * {sum_j} = {numerator:.6f}")
    print(f"The denominator of the expression for <k> = 1 + k*x = 1 + {k}*{x:.6f} = {denominator:.6f}")
    print(f"\nAverage number of adsorbed layers per site = {numerator:.6f} / {denominator:.6f} = {average_layers:.6f}")
    print("<<<" + str(average_layers) + ">>>")

solve_adsorption()