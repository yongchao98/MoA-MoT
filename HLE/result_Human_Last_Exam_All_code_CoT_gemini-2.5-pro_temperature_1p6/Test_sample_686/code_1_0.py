def solve_ising_susceptibility():
    """
    This function prints the derivation and the final formula for the magnetic susceptibility.
    """
    # Define variables as strings for a clear representation of the formula.
    chi = "\\chi"
    beta = "\\beta"
    c = "c"
    l = "l"
    C_l = "C_l"
    m_0 = "m_0"
    B_l = "B_l"
    T_prop = "T"
    N = "N"

    # Step 1: Start with the given formula for susceptibility
    print("Step 1: The given formula for the magnetic susceptibility is:")
    print(f"{chi} = {beta} * \\sum_{{{l}=1}}^{{\\infty}} {c}({c}-1)^{{{l}-1}} {C_l}\n")

    # Step 2: Substitute the definition of the connected correlation
    print("Step 2: We substitute the relation for the connected correlation:")
    print(f"{C_l} = (1/{beta}) * d<{m_0}>/d{B_l}")
    print("This gives:")
    print(f"{chi} = \\sum_{{{l}=1}}^{{\\infty}} {c}({c}-1)^{{{l}-1}} d{m_0}/d{B_l}\n")

    # Step 3 & 4: Evaluate the derivative by considering perturbation propagation
    print("Step 3 & 4: We evaluate the derivative d{m_0}/d{B_l}.")
    print("A perturbation dB_l at site l propagates along the path to site 0.")
    print("The effect on the magnetization at site 0 is found using the chain rule.")
    print(f"The influence decays at each step by a factor {T_prop}, the propagation transfer factor.")
    print(f"This factor T is given by T = \\tanh(\\beta J) * (1 - m_cav^2) / (1 - \\tanh^2(\\beta J)m_cav^2),")
    print(f"where m_cav is the cavity magnetization.\n")

    # Step 5: Express C_l
    print("Step 5: Following the propagation over a path of length l, we find:")
    print(f"d{m_0}/d{B_l} = {beta} * (1 - {m_0}^2) * {T_prop}^{l}")
    print("Therefore, the connected correlation is:")
    print(f"{C_l} = (1 - {m_0}^2) * {T_prop}^{l}\n")

    # Step 6: Substitute C_l back and sum the series
    print("Step 6: Substituting C_l into the expression for chi and summing the geometric series:")
    print(f"{chi} = \\sum_{{{l}=1}}^{{\\infty}} {c}({c}-1)^{{{l}-1}} (1 - {m_0}^2) * {T_prop}^{l}")
    print(f"{chi} = {c} * (1 - {m_0}^2) * {T_prop} * \\sum_{{{l}=1}}^{{\\infty}} (({c}-1){T_prop})^{{l-1}}")
    print("The sum of the geometric series is 1 / (1 - (c-1)T). This gives:")
    print(f"{chi} = ({beta} * {c} * (1 - {m_0}^2) * {T_prop}) / (1 - ({c}-1)*{T_prop})\n")

    # Step 7: Use the provided constant N for simplification
    print("Step 7: We use the provided constant N for brevity.")
    print(f"Given {N} = {beta} * ({c}(1-{m_0}^2)) / ({c}-1), which means {beta}*c*(1-{m_0}^2) = N*({c}-1).")
    print("Substituting this into the expression for chi yields the final answer:\n")
    
    # Final Answer
    final_chi_expression = f"{chi} = ( {N} * ({c}-1) * {T_prop} ) / ( {beta} * (1 - ({c}-1) * {T_prop}) )"
    print("Final Answer:")
    print(final_chi_expression)


solve_ising_susceptibility()