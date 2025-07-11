def display_green_function_dependence():
    """
    This function demonstrates the functional dependence of the bare Green's function G_0
    on the single-particle energy eigenvalues ϵ_k by printing its formula.
    """

    # Define the symbols used in the equation as strings
    g0_symbol = "G_0(k, ω)"
    equals = "="
    numerator = "1"
    omega = "ω"
    epsilon_k = "ϵ_k"
    infinitesimal = "iδ"

    # Explain the components of the equation
    print("The bare Green's function, G_0, has an inverse dependence on the single-particle energy, ϵ_k.")
    print("The formula is constructed from the following parts:")
    print(f"\n1. The function itself: {g0_symbol}")
    print(f"2. A numerator, which is: {numerator}")
    print(f"3. A denominator composed of three terms:")
    print(f"   - Frequency: {omega}")
    print(f"   - Single-particle energy: {epsilon_k}")
    print(f"   - An infinitesimal for causality: {infinitesimal}")

    # Assemble and print the final equation
    # The final equation shows G_0 = 1 / (ω - ϵ_k + iδ)
    final_equation = f"{g0_symbol} {equals} {numerator} / ({omega} - {epsilon_k} + {infinitesimal})"

    print("\n--- Final Functional Dependence ---")
    print(final_equation)
    print("---------------------------------")
    print("\nThis equation shows that G_0(k, ω) has a pole when the frequency ω is equal to the particle's energy ϵ_k.")

# Execute the function to display the answer
display_green_function_dependence()