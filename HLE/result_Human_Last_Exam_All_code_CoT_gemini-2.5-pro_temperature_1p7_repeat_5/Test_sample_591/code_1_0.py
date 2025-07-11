import sys

def solve_task():
    """
    This script provides the definition of the parameter κ from the dendritic plasticity model.
    It includes a function to calculate κ and demonstrates its use with example values.
    """

    def calculate_kappa(phi, mu, rho, tau_u, tau_v):
        """
        Calculates the parameter kappa based on the model's parameters.

        Args:
            phi (float): Presynaptic scaling constant.
            mu (float): Mean firing rate.
            rho (float): Offset constant in the Hebbian rule.
            tau_u (float): Time constant for the postsynaptic accumulator.
            tau_v (float): Time constant for the presynaptic accumulator.

        Returns:
            float: The calculated value of kappa.
            Returns an error string if input parameters are invalid.
        """
        # Ensure parameters that will be in the denominator are not zero
        if phi == 0 or tau_u == 0 or tau_v == 0:
            return "Error: parameters phi, tau_u, and tau_v cannot be zero."

        # Calculate kappa using the derived formula
        term1 = mu + rho / phi
        term2 = 1/tau_u + 1/tau_v
        kappa = -term1 * term2
        return kappa

    # Explain the definition of kappa
    print("Based on the stability analysis of the synaptic weight dynamics, the parameter κ in the expression for critical correlation is defined by the following formula:")
    print()
    print("κ = -(μ + ρ/φ) * (1/τ_u + 1/τ_v)")
    print()
    print("where:")
    print("  μ: mean firing rate of neurons")
    print("  ρ: offset constant for plasticity")
    print("  φ: presynaptic scaling constant")
    print("  τ_u: postsynaptic accumulator time constant")
    print("  τ_v: presynaptic accumulator time constant")
    print("-" * 30)

    # Demonstrate the calculation with example numerical values
    # These values are for demonstration purposes only.
    phi_val = 1.0
    mu_val = 10.0   # in Hz
    rho_val = -15.0 # dimensionless
    tau_u_val = 20.0  # in ms
    tau_v_val = 15.0  # in ms

    kappa_value = calculate_kappa(phi_val, mu_val, rho_val, tau_u_val, tau_v_val)

    print("Example Calculation:")
    print(f"Given the parameters:")
    print(f"  μ = {mu_val}")
    print(f"  ρ = {rho_val}")
    print(f"  φ = {phi_val}")
    print(f"  τ_u = {tau_u_val}")
    print(f"  τ_v = {tau_v_val}")
    print()

    # Show the calculation step-by-step with the example numbers
    term1_val = mu_val + rho_val / phi_val
    term2_val = 1/tau_u_val + 1/tau_v_val

    print("The final equation for κ with these numbers is:")
    print(f"κ = -({mu_val} + {rho_val} / {phi_val}) * (1 / {tau_u_val} + 1 / {tau_v_val})")
    print(f"κ = -({term1_val}) * ({term2_val:.4f})")
    print(f"κ = {kappa_value:.4f}")

solve_task()