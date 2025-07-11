def solve_nsvz_fixed_point(Nc, Nf):
    """
    Calculates the anomalous dimension at the infrared fixed point
    of N=1 SQCD using the NSVZ beta function relation.
    The fixed point exists when the numerator of the beta function is zero.

    Args:
        Nc (int): The number of colors (e.g., 3 for SU(3)).
        Nf (int): The number of quark flavors.
    """
    # The NSVZ beta function for N=1 SQCD with gauge group SU(Nc) and Nf flavors
    # (in the fundamental and anti-fundamental representations) is given by:
    # beta(g) = -g^3/(16*pi^2) * [3*Nc - Nf*(1-gamma(g))] / [1 - g^2*Nc/(8*pi^2)]
    # A non-trivial infrared fixed point g* exists if beta(g*) = 0.
    # This requires the numerator of the fraction to be zero.
    
    # We first check for asymptotic freedom, a necessary condition for the theory
    # to flow from weak coupling in the UV to a potentially strong-coupling
    # fixed point in the IR.
    one_loop_coeff = 3 * Nc - Nf
    if one_loop_coeff <= 0:
        print(f"For Nc={Nc} and Nf={Nf}, the theory is not asymptotically free.")
        print("It does not flow to an interacting infrared fixed point from the UV.")
        return

    # At the fixed point, the numerator is zero: 3*Nc - Nf*(1 - gamma) = 0
    # We can solve this equation for the anomalous dimension 'gamma' of the matter superfield.
    # Note: The result depends on the renormalization scheme. The NSVZ formula as
    # written is valid in the "NSVZ scheme". In other schemes, like DRED, the
    # relation between beta and gamma is modified.
    
    # Solving for gamma:
    # Nf*(1 - gamma) = 3*Nc
    # 1 - gamma = 3*Nc / Nf
    # gamma = 1 - (3*Nc / Nf)
    gamma_at_fixed_point = 1.0 - (3.0 * Nc / Nf)

    print("The condition for an infrared fixed point from the NSVZ beta function is that its numerator vanishes.")
    print(f"For N=1 SQCD with SU({Nc}) gauge group and {Nf} flavors, the equation is:")
    print("3*Nc - Nf*(1 - gamma) = 0")
    print("\nSolving this equation for the anomalous dimension 'gamma':")
    
    # Print each number in the final equation step-by-step
    print(f"3 * {Nc} - {Nf} * (1 - gamma) = 0")
    print(f"{3*Nc} - {Nf} + {Nf}*gamma = 0")
    print(f"{Nf}*gamma = {Nf - 3*Nc}")
    print(f"gamma = ({Nf} - {3*Nc}) / {Nf}")
    print(f"gamma = {1 - (3*Nc)}/ {Nf}")
    print(f"gamma = 1 - {3*Nc / Nf:.4f}")
    print(f"gamma = {gamma_at_fixed_point:.4f}")

# Let's use an example from the conformal window of N=1 SQCD.
# For example, Nc=3, Nf=5.
print("--- Calculating fixed point for Nc=3, Nf=5 ---")
solve_nsvz_fixed_point(Nc=3, Nf=5)