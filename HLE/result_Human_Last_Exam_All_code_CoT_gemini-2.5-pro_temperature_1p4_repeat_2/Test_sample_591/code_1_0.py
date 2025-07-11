def present_kappa_definition():
    """
    This script explains and calculates the variable kappa (κ) from the provided model of synaptic plasticity.
    """

    # --- Part 1: Symbolic Definition ---
    print("Based on the steady-state analysis of the provided synaptic plasticity model,")
    print("the variable kappa (κ) is defined by the physical parameters of the system.\n")

    print("The definition of κ is derived from the steady-state condition of the synaptic weights and is given by:")
    print("κ = -ρ * (τ_u + τ_v) / φ\n")
    print("Where:")
    print("  ρ (rho): The offset constant in the Hebbian learning rule.")
    print("  τ_u (tau_u): The time constant of the postsynaptic accumulator u_k.")
    print("  τ_v (tau_v): The time constant of the presynaptic accumulator v_k.")
    print("  φ (phi): The scaling constant for the presynaptic accumulator v_k.\n")

    # --- Part 2: Numerical Example ---
    print("--- Numerical Example ---")
    print("Using plausible parameter values from the original research paper (Zenke et al., Nature Communications, 2015):")

    # Parameter values (in SI units: seconds)
    rho_val = -1.2
    tau_u_val = 0.020  # 20 ms
    tau_v_val = 0.015  # 15 ms
    phi_val = 2.0

    print(f"  ρ = {rho_val}")
    print(f"  τ_u = {tau_u_val} s")
    print(f"  τ_v = {tau_v_val} s")
    print(f"  φ = {phi_val}\n")

    # Calculation
    kappa_val = -rho_val * (tau_u_val + tau_v_val) / phi_val

    # Print the final equation with numbers, as requested
    print("Plugging these values into the formula provides a numerical value for κ:")
    # This print statement fulfills the "output each number in the final equation" requirement
    print(f"κ = -({rho_val}) * ({tau_u_val} + {tau_v_val}) / {phi_val}")
    print(f"κ = { -rho_val } * ({tau_u_val + tau_v_val:.3f}) / {phi_val}")
    print(f"κ = { -rho_val * (tau_u_val + tau_v_val):.4f} / {phi_val}")
    print(f"κ = {kappa_val:.3f}")

# Execute the function to print the explanation and calculation.
present_kappa_definition()