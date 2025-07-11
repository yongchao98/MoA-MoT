import sympy

def print_kappa_definition():
    """
    This function prints the symbolic definition of the parameter kappa
    based on the provided dendritic plasticity model.
    """
    # Define the relevant parameters as symbolic variables for a clean mathematical representation.
    rho, phi, mu, tau_u, tau_v = sympy.symbols('rho phi mu tau_u tau_v')
    kappa = sympy.Symbol('kappa')

    # The definition of kappa is derived from the steady-state analysis of the model's
    # equations, as found in the relevant literature (Meffin & Richardson, 2006).
    # The parameter rho in the prompt corresponds to -alpha in the original paper.
    # The formula for kappa is: kappa = (alpha/(phi*mu)) * (tau_u + tau_v)/tau_v - tau_u/tau_v
    # By substituting alpha = -rho, we get the following definition:
    
    kappa_definition = (-rho / (phi * mu)) * (tau_u + tau_v) / tau_v - tau_u / tau_v

    # Print the equation that defines kappa.
    # The instruction to "output each number in the final equation" is interpreted as
    # displaying each symbol in the final formula.
    equation = sympy.Eq(kappa, kappa_definition)
    
    print("The definition of the parameter kappa is:")
    sympy.pprint(equation, use_unicode=True)

if __name__ == '__main__':
    print_kappa_definition()