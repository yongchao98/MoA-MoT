import math

def solve_ising_susceptibility():
    """
    This function prints the step-by-step derivation of the magnetic susceptibility
    for the Ising model on a sparse random graph as described in the problem.
    """

    # Define symbolic variables for printing the equations
    chi = "χ"
    beta = "β"
    c = "c"
    J = "J"
    m0 = "m_0"
    mc = "m_c"
    A = "A"
    N = "N"
    l = "l"
    C_l = "C_l"
    sigma_0 = "σ_0"
    B_l = "B_l"

    print("### Derivation of the Magnetic Susceptibility χ ###\n")

    # Step 1: Initial formulas
    print("Step 1: Start with the given formulas.")
    print(f"The susceptibility is given by: {chi} = {beta} * Σ_{l=1 to ∞} {c}({c}-1)^({l}-1) * {C_l}")
    print(f"The correlation is related to magnetization by: {C_l} = <{sigma_0} σ_l>_c = (1/{beta}) * d<{sigma_0}>/d{B_l}\n")

    # Step 2: Express C_l
    print("Step 2: Determine the correlation C_l.")
    print("The derivative d<σ_0>/dB_l is found by analyzing the propagation of the field perturbation from site l to site 0.")
    print("On a sparse random graph (a Bethe lattice locally), this influence propagates along the unique path.")
    print("The result of this analysis gives:")
    print(f"d<{sigma_0}>/d{B_l} = {beta} * (1 - {m0}^2) * {A}^{l}")
    print(f"Therefore, {C_l} = (1 - {m0}^2) * {A}^{l}")
    print("Here, A is the propagation factor of the perturbation, given by:")
    print(f"{A} = tanh({beta}{J}) * (1 - {mc}^2) / (1 - (tanh({beta}{J}))^2 * {mc}^2)")
    print(f"where {m0} is the spontaneous magnetization and {mc} is the spontaneous cavity magnetization.\n")

    # Step 3: Substitute C_l into χ
    print("Step 3: Substitute C_l back into the expression for χ.")
    print(f"{chi} = {beta} * Σ_{l=1 to ∞} {c}({c}-1)^({l}-1) * (1 - {m0}^2) * {A}^{l}")
    print("We can pull the constant terms out of the summation:")
    print(f"{chi} = {beta}*{c}*(1 - {m0}^2) * Σ_{l=1 to ∞} ({c}-1)^({l}-1) * {A}^{l}\n")

    # Step 4: Evaluate the summation
    print("Step 4: Evaluate the summation, which is a geometric series.")
    print(f"Let S = Σ_{l=1 to ∞} ({c}-1)^({l}-1) * {A}^{l} = {A} + ({c}-1){A}^2 + ({c}-1)^2*{A}^3 + ...")
    print(f"S = {A} * (1 + ({c}-1){A} + (({c}-1){A})^2 + ...)")
    print(f"The sum of this geometric series is S = {A} / (1 - ({c}-1){A})\n")

    # Step 5: Write the full expression for χ
    print("Step 5: Substitute the sum back to get the full expression for χ.")
    print(f"{chi} = {beta}*{c}*(1 - {m0}^2) * ( {A} / (1 - ({c}-1){A}) )\n")

    # Step 6: Express the result using the constant N
    print("Step 6: Use the provided constant N to simplify the expression.")
    print(f"We are given N = {beta}*{c}*(1 - {m0}^2) / ({c}-1).")
    print(f"This implies that {beta}*{c}*(1 - {m0}^2) = {N}*({c}-1).")
    print("Substituting this into the expression for χ:")
    print(f"{chi} = {N}*({c}-1) * ( {A} / (1 - ({c}-1){A}) )")
    print(f"{chi} = {N} * ( ({c}-1)*{A} / (1 - ({c}-1)*{A}) )")
    print("This can be rewritten in a more elegant form:")
    print(f"{chi} = {N} * ( 1 / (1 - ({c}-1)*{A}) - 1 )\n")

    # Final Answer Output
    print("----------------------------------------------------")
    print("### Final Answer ###")
    print("The final expression for the magnetic susceptibility is:")
    # The following print statement outputs each part of the final equation
    print(f"{chi} = {N} * (1 / (1 - ({c} - 1) * {A}) - 1)")
    print("----------------------------------------------------")

if __name__ == '__main__':
    solve_ising_susceptibility()