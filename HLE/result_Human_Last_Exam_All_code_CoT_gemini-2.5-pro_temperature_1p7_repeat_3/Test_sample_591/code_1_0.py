import sys

# The user is asking for the definition of the parameter 'kappa'
# from a model of dendritic plasticity. The definition is a formula
# composed of other parameters from the model.

# The relevant parameters mentioned or implied in the text are:
# phi: The scaling constant for the presynaptic accumulator.
# rho: The offset constant in the Hebbian learning rule.
# tau_w: The time constant for the synaptic efficacy change.
# tau_u: The time constant for the postsynaptic accumulator.
# mu: The mean firing rate, assumed equal for all synapses at steady state.
# W_thr: The weight threshold for synaptic turnover.

# The script below will print the definition of kappa using these parameters.

def print_kappa_definition():
    """
    Prints the formula defining kappa and describes its component variables.
    """
    print("The definition of the term kappa (κ) in the expression for critical correlation is:")
    print("")

    # The final equation defining kappa
    final_equation = "kappa = -(phi / rho) * (tau_w * mu / tau_u) * W_thr"
    
    # We must output each component of the final equation as per the instructions.
    # For a definitional formula, this means printing the formula itself.
    print(final_equation)
    print("")

    print("Where each variable in the formula represents a parameter from the model:")
    print("  phi   : The scaling constant from the presynaptic accumulator equation.")
    print("  rho   : The offset constant from the Hebbian learning rule.")
    print("  tau_w : The time constant for the change in synaptic efficacy.")
    print("  mu    : The mean firing rate of synapses.")
    print("  tau_u : The time constant for the postsynaptic accumulator.")
    print("  W_thr : The fixed threshold for synaptic efficacy below which a synapse is removed.")

# Execute the function to provide the definition.
print_kappa_definition()
