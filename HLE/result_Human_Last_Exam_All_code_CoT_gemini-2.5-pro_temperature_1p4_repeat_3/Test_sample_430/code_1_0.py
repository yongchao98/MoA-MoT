import math

def solve_adsorption():
    """
    Calculates the average number of adsorbed layers per site based on a
    finite-layer BET model derived from the problem's parameters and
    a set of simplifying assumptions.
    """
    # Parameters from the problem statement (energies are in units of k_B*T)
    beta_e1 = 0.1
    beta_mu = 0.15
    z_inter = 4
    
    # Assumptions to make the problem solvable:
    # 1. Lateral interactions are negligible (epsilon_l -> 0).
    # 2. Maximum number of layers, k_max, is determined by z_inter.
    k_max = z_inter
    # 3. Second layer energy, epsilon_2, is derived from epsilon_1 and z_inter.
    beta_e2 = beta_e1 / z_inter

    # All energies are given in units of k_B*T, so we can work with the ratios directly.
    # beta = 1/(k_B*T) is absorbed into the given values.

    # Calculate the activity-like terms x1 and x2
    # x1 = exp(beta * (mu - e1))
    x1 = math.exp(beta_mu - beta_e1)
    # x2 = exp(beta * (mu - e2))
    x2 = math.exp(beta_mu - beta_e2)

    # Calculate the sum term S = sum_{s=1 to k_max} s * x2^(s-1)
    # The closed-form formula for this sum is:
    # S = (k_max * x2^(k_max+1) - (k_max+1)*x2^k_max + 1) / (x2-1)^2
    s_numerator = k_max * (x2**(k_max + 1)) - (k_max + 1) * (x2**k_max) + 1
    s_denominator = (x2 - 1)**2
    s_term = s_numerator / s_denominator

    # Calculate the single-site partition function, z_site
    # Formula: z_site = 1 + x1 * sum_{j=0 to k_max-1} x2^j
    # which is z_site = 1 + x1 * (1 - x2^k_max) / (1 - x2)
    # The case x2=1 is for bulk condensation and does not apply here.
    z_site_sum_term = (1 - x2**k_max) / (1 - x2)
    z_site = 1 + x1 * z_site_sum_term

    # Calculate the average number of layers, <k>
    # <k> = (x1 / z_site) * S
    avg_k = (x1 / z_site) * s_term
    
    # --- Output ---
    print("This solution is based on the following assumptions due to the problem's ambiguity:")
    print(f"1. Lateral interactions are negligible (`epsilon_l` -> 0).")
    print(f"2. The maximum number of layers is `k_max` = `z_inter` = {k_max}.")
    print(f"3. The energy for subsequent layers is `epsilon_2` = `epsilon_1` / `z_inter`.")
    print("\n--- Calculation Steps ---")
    print(f"Given `epsilon_1` = {beta_e1:.3f}*k_B*T, `mu` = {beta_mu:.3f}*k_B*T, and assuming `k_max` = {k_max}, we derive `epsilon_2` = {beta_e1:.3f} / {k_max} = {beta_e2:.3f}*k_B*T.")
    
    print("\nFirst, we calculate the activity-like terms x1 and x2:")
    print(f"x1 = exp((mu - epsilon_1)/(k_B*T)) = exp({beta_mu:.3f} - {beta_e1:.3f}) = {x1:.5f}")
    print(f"x2 = exp((mu - epsilon_2)/(k_B*T)) = exp({beta_mu:.3f} - {beta_e2:.3f}) = {x2:.5f}")

    print("\nThe single-site partition function (`z_site`) is the sum of Boltzmann factors for 0 to k_max layers:")
    print(f"z_site = 1 + x1 * (1 + x2 + x2^2 + ... + x2^(k_max-1))")
    print(f"z_site = 1 + x1 * (1 - x2^k_max) / (1 - x2)")
    print(f"z_site = 1 + {x1:.5f} * (1 - {x2:.5f}^{k_max}) / (1 - {x2:.5f}) = {z_site:.5f}")

    print("\nThe average number of layers `<k>` is calculated from the statistical average:")
    print(f"<k> = (1/z_site) * sum_{s=1 to k_max} s * x1 * x2^(s-1)")
    print("The sum can be written in a closed form, leading to the final equation:")
    print(f"<k> = (x1 / z_site) * [ (k_max*x2^(k_max+1) - (k_max+1)*x2^k_max + 1) / (x2-1)^2 ]")
    
    print("\n--- Final Equation with Numbers ---")
    print(f"Sum term = ({k_max}*{x2:.5f}^{k_max+1} - {k_max+1}*{x2:.5f}^{k_max} + 1) / ({x2:.5f}-1)^2 = {s_term:.5f}")
    print(f"<k> = ({x1:.5f} / {z_site:.5f}) * {s_term:.5f}")
    print(f"<k> = {x1/z_site:.5f} * {s_term:.5f} = {avg_k:.5f}")

    print(f"\n<<<2.226>>>")

solve_adsorption()