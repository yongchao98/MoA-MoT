import math

def calculate_average_layers():
    """
    Calculates the average number of adsorbed layers per site based on a
    stabilized multi-layer adsorption model.

    The core challenge is the positive chemical potential (mu > 0), which
    in standard models leads to infinite adsorption. To find a finite answer,
    we must assume a model that ensures stability (activity X < 1). This forces
    the interaction energy for subsequent layers (eps_L) to be repulsive.

    The plan is as follows:
    1. Define parameters in units of k_B*T from the problem statement.
    2. Hypothesize a repulsive energy eps_L = -z_inter * eps_1 to ensure stability.
    3. Assume lateral interactions are negligible as k->infinity.
    4. Use the BET formula for infinite layers in the stable regime to calculate <k>.
    5. Print the equation with all numerical values and the final result.
    """
    # 1. Define parameters from the problem (in units of k_B*T)
    mu_kT = 0.15
    eps1_kT = 0.1
    z_inter = 4.0

    # 2. Define the hypothesized repulsive energy for subsequent layers
    # This is the key assumption to make the problem solvable.
    eps_L_kT = -z_inter * eps1_kT

    # 3. Calculate the activity terms C and X
    # C is for the first layer, X is for subsequent layers.
    C = math.exp(mu_kT + eps1_kT)
    X = math.exp(mu_kT + eps_L_kT)

    # 4. Calculate the average number of layers <k> using the BET formula
    # for an infinite number of layers where X < 1.
    denominator = (1 - X) * (1 - X + C)
    average_k = C / denominator

    # 5. Print the breakdown of the calculation as requested.
    print(f"Based on the interpretation to ensure a stable system:")
    print(f"Assumed energy for layers > 1: eps_L/kT = -z_inter * eps1/kT = {eps_L_kT:.2f}")
    print(f"Activity for the first layer: C = exp({mu_kT} + {eps1_kT}) = {C:.4f}")
    print(f"Activity for subsequent layers: X = exp({mu_kT} + {eps_L_kT:.2f}) = {X:.4f}")
    
    print(f"\nThe final equation with numerical values is:")
    print(f"<k> = {C:.4f} / ((1 - {X:.4f}) * (1 - {X:.4f} + {C:.4f}))")

    print(f"\nThe calculated average number of adsorbed layers per site is:")
    print(f"<k> = {average_k:.2f}")

calculate_average_layers()