import math

def display_kappa_definition():
    """
    This function prints the definition of the parameter kappa (κ) from the critical
    correlation formula, based on a derivation from the provided plasticity model.
    """

    # The parameter kappa represents a dimensionless combination of the model's
    # physical parameters. It essentially sets the critical correlation required
    # for a synapse to stabilize when it is part of a large group of synapses.
    # It can be understood as a ratio between mean-driven and fluctuation-driven
    # components of the learning rule.

    # The final equation for kappa is derived as κ = -(μ*(φ*μ + ρ))/(φ*σ^2).
    # The code below will print this formula and its components.
    
    # Define the symbols as strings for printing
    kappa = "κ"
    mu = "μ"
    phi = "φ"
    rho = "ρ"
    sigma_sq = "σ^2"

    print("The definition of the parameter κ in the given expression is:")
    
    # Printing the final equation by describing its components
    # to fulfill the requirement of outputting each part of the equation.
    numerator_part_1 = mu
    numerator_part_2 = f"({phi}*{mu} + {rho})"
    denominator_part_1 = phi
    denominator_part_2 = sigma_sq
    
    print(f"{kappa} = - ( {numerator_part_1} * {numerator_part_2} ) / ( {denominator_part_1} * {denominator_part_2} )")

    print("\nWhere the symbols represent the following model parameters:")
    print(f"{mu} (mu): The mean activity (e.g., firing rate) of presynaptic neurons.")
    print(f"{sigma_sq} (sigma^2): The variance of the presynaptic neurons' activity.")
    print(f"{phi} (phi): The scaling constant for the presynaptic accumulator.")
    print(f"{rho} (rho): The offset constant in the Hebbian learning rule.")

if __name__ == '__main__':
    display_kappa_definition()
