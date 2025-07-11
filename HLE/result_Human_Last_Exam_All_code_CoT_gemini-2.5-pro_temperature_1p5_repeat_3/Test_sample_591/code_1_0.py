import sys

def get_kappa_definition():
    """
    This function prints the definition of kappa (κ) based on the provided model
    of local dendritic plasticity. The definition is derived from the source literature
    for this specific model, as it is not explicitly stated in the user's prompt.
    """
    # Define the symbols for clarity using their text representations.
    # These parameters are described in the model.
    rho = "ρ"         # offset constant
    tau_w = "τ_w"     # time constant for synaptic efficacy
    W_thr = "W_thr"   # synaptic efficacy threshold for removal
    phi = "φ"         # presynaptic scaling constant
    mu = "μ"          # mean firing rate of synapses
    tau_u = "τ_u"     # time constant for postsynaptic accumulator
    tau_v = "τ_v"     # time constant for presynaptic accumulator

    # Construct the numerator and denominator strings for the kappa equation
    # The instruction "output each number in the final equation" seems to be
    # a general template instruction that doesn't apply here, as there are no
    # numerical values provided. The code will instead output the symbolic equation.
    kappa_numerator_str = f"-({rho} * {tau_w} * {W_thr})"
    kappa_denominator_str = f"{phi} * {mu} * {tau_u} * {tau_v}"
    
    # Python 3.7+ is needed for the unicode characters.
    # Fallback for older versions.
    if sys.version_info < (3, 7):
        kappa_numerator_str = f"-(rho * tau_w * W_thr)"
        kappa_denominator_str = f"phi * mu * tau_u * tau_v"

    print("The definition of κ is a dimensionless parameter that groups several constants from the model:")
    print("\n")
    # Pretty print the equation
    print(f"       {kappa_numerator_str}")
    print(f"  κ =  {'-' * (len(kappa_numerator_str) + 2)}")
    print(f"       {kappa_denominator_str}")
    print("\n")
    print("Where the parameters are:")
    print(f"  {rho} (rho): The offset constant in the Hebbian learning rule.")
    print(f"  {tau_w} (tau_w): The time constant for the synaptic efficacy (w).")
    print(f"  {W_thr} (W_thr): The efficacy threshold for synapse removal.")
    print(f"  {phi} (phi): The scaling constant for the presynaptic accumulator (v).")
    print(f"  {mu} (mu): The mean firing rate of the synapses.")
    print(f"  {tau_u} (tau_u): The time constant for the postsynaptic accumulator (u).")
    print(f"  {tau_v} (tau_v): The time constant for the presynaptic accumulator (v).")

# Execute the function to print the definition
get_kappa_definition()
