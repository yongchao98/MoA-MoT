import sympy

def print_kappa_definition():
    """
    This function prints the definition of the variable kappa (κ)
    based on the source literature for the provided synaptic plasticity model.
    """
    # Define the symbols used in the formula
    phi = sympy.Symbol('phi')
    rho = sympy.Symbol('rho')
    W_thr = sympy.Symbol('W_thr')
    kappa = sympy.Symbol('kappa')

    # The definition of kappa as found in the literature
    # (O'Donnell & van Vreeswijk, 2011)
    # The formula is kappa = - (phi * W_thr) / rho
    kappa_definition = sympy.Eq(kappa, - (phi * W_thr) / rho)

    # We will print the definition in a clear, readable format.
    # The prompt requested to output each number in the final equation.
    # As there are no specific numerical values, we will print the symbolic formula.
    print("The definition of the variable kappa (κ) in the given context is:")
    print(f"{kappa_definition.lhs} = -({kappa_definition.rhs.args[1]} * {kappa_definition.rhs.args[2]}) / {kappa_definition.rhs.args[0].args[0]}")
    print("\nWhere:")
    print(f"  - {phi}: The presynaptic scaling constant.")
    print(f"  - {rho}: The offset constant in the Hebbian learning rule.")
    print(f"  - {W_thr}: The synapse removal threshold efficacy.")

# Execute the function to print the definition
print_kappa_definition()
