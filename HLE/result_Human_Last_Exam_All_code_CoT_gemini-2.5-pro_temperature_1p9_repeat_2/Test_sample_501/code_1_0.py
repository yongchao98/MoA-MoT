import math

def calculate_isolated_polymer_force(n, l, E0, x):
    """
    Calculates the restorative force of a thermally isolated freely jointed chain.

    Args:
        n (int): Number of links in the polymer chain.
        l (float): Length of each link.
        E0 (float): Kinetic energy of the polymer at zero extension (in Joules).
        x (float): The end-to-end extension of the polymer.
    """

    print(f"Calculating the force for a polymer with the following parameters:")
    print(f"  Number of links (n) = {n}")
    print(f"  Length of each link (l) = {l}")
    print(f"  Initial kinetic energy (E(0)) = {E0:.2e} J")
    print(f"  Extension (x) = {x}\n")

    # The derived force law is f = (2 * E(0) * x / (n^2 * l^2)) * exp(x^2 / (n^2 * l^2))
    # We will print the equation with the numbers substituted.
    
    # Calculate terms for the equation string
    n_squared_l_squared = n**2 * l**2
    x_squared = x**2
    
    # Construct the equation string with substituted values
    equation_string_numerator_1 = f"2 * {E0:.2e} * {x}"
    equation_string_denominator_1 = f"{n}**2 * {l}**2"
    
    equation_string_exp_numerator = f"{x}**2"
    equation_string_exp_denominator = f"{n}**2 * {l}**2"

    print("The force law for a thermally isolated polymer is:")
    print("f = (2 * E(0) * x / (n^2 * l^2)) * exp(x^2 / (n^2 * l^2))\n")
    
    print("Substituting the values into the equation:")
    # Print the detailed calculation steps as an equation
    print(f"f = ( {equation_string_numerator_1} / ({equation_string_denominator_1}) ) * exp( {equation_string_exp_numerator} / ({equation_string_exp_denominator}) )")

    # Calculate the final force value
    pre_factor = (2 * E0 * x) / n_squared_l_squared
    exp_term = math.exp(x_squared / n_squared_l_squared)
    force = pre_factor * exp_term

    # Print the intermediate values to show the full calculation
    print(f"f = ( {pre_factor:.3e} ) * exp( {x_squared / n_squared_l_squared:.3e} )")
    print(f"f = ( {pre_factor:.3e} ) * ( {exp_term:.5f} )")
    
    print("\nFinal calculated force:")
    print(f"f = {force:.3e} N")


# Example parameters for a polymer chain
# n should be large for the Gaussian approximation to hold
num_links = 1000
# l can be in any unit, but we'll use meters for SI units
link_length = 1.0e-9  # 1 nm
# x should be small compared to n*l
extension = 30e-9     # 30 nm, a reasonable extension
# E(0) can be estimated from an initial temperature T, e.g., 300 K (room temp)
# E(0) = (3/2) * n * k_B * T, where k_B is the Boltzmann constant
boltzmann_constant = 1.380649e-23  # J/K
initial_temp = 300  # K
initial_energy = (3/2) * num_links * boltzmann_constant * initial_temp

calculate_isolated_polymer_force(n=num_links, l=link_length, E0=initial_energy, x=extension)
