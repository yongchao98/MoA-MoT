import sys

def solve_susceptibility():
    """
    This function programmatically prints the derivation and final formula for the
    magnetic susceptibility chi of an Ising model on a sparse random graph.
    """
    # Define symbolic variables as strings for clarity in the output
    beta = "beta"
    c = "c"
    J = "J"
    B = "B"
    l = "l"
    m0_sq = "m_0^2"
    mcav = "m_cav"
    mcav_sq = f"{mcav}^2"
    one = "1"
    
    # Header
    print("Derivation of the magnetic susceptibility chi:")
    print("=" * 45)

    # Step 1: Define the propagation factor T
    # T is the factor by which the effect of a field perturbation is attenuated
    # as it passes from one site to the next along a chain. It is derived from
    # the message-passing equations under homogeneous conditions.
    print("Step 1: Define the propagation factor T")
    print("Let T be the factor describing the propagation of a perturbation along the spin chain.")
    T_numerator = f"tanh({beta}*{J}) * ({one} - {mcav_sq})"
    T_denominator = f"{one} - (tanh({beta}*{J}))^2 * {mcav_sq}"
    print(f"T = ({T_numerator}) / ({T_denominator})")
    print(f"where {mcav} is the cavity magnetization satisfying the system's self-consistency equations.")
    print("-" * 45)

    # Step 2: Express the correlation C_l
    # By applying the chain rule to d<sigma_0>/dB_l, we find the correlation C_l.
    print("Step 2: Calculate the correlation C_l")
    print("The correlation C_l = <sigma_0 sigma_l>_c is found by propagating the perturbation from site l to 0.")
    C_l_expression = f"({one} - {m0_sq}) * T^{l}"
    print(f"C_l = (1/beta) * d<sigma_0>/dB_l = {C_l_expression}")
    print(f"where m_0 is the magnetization at site 0.")
    print("-" * 45)

    # Step 3: Evaluate the sum for chi
    print("Step 3: Evaluate the sum for susceptibility")
    print("Substituting C_l into the formula for chi and evaluating the geometric series:")
    print(f"chi = {beta} * Sum_{{l=1..inf}} [ c*({c}-{one})^(l-{one}) * ({one}-{m0_sq})*T^{l} ]")
    print("This simplifies to a closed-form expression.")
    print("-" * 45)
    
    # Step 4: Final expression for chi using the given constant N
    print("Step 4: Final simplified expression for chi")
    print("Using the definition from the hint, N = beta * c * (1 - m_0^2) / (c-1), the expression becomes:")
    
    # Define the components of the final formula
    N_var = "N"
    c_minus_1 = f"({c} - {one})"
    T_var = "T"
    
    numerator = f"{N_var} * {c_minus_1} * {T_var}"
    denominator = f"{one} - {c_minus_1} * {T_var}"
    
    print(f"chi = ({numerator}) / ({denominator})")
    print("\n")
    
    # Fulfilling the requirement to "output each number in the final equation"
    print("Here are the individual components and numbers of the final equation:")
    print(f"Term 'N': This is the constant N = beta * c * (1 - m_0^2) / (c-1).")
    print(f"Term '{c_minus_1}': The effective branching factor, connectivity 'c' minus the number '{one}'.")
    print(f"Term 'T': The propagation factor T.")
    print(f"The number in the denominator is '{one}'.")


if __name__ == '__main__':
    solve_susceptibility()
