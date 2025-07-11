import math

def calculate_isolated_polymer_force(x, l, n, E0):
    """
    Calculates the restorative force of a thermally isolated, freely jointed polymer chain.

    This function is based on the microcanonical ensemble where total entropy is conserved.

    Args:
        x (float): The end-to-end separation of the polymer.
        l (float): The length of a single monomer link.
        n (int): The number of links in the polymer chain.
        E0 (float): The kinetic energy of the polymer at zero extension (x=0).

    Returns:
        float: The force exerted by the polymer. A negative value indicates a
               restoring force, pulling the ends together.
    """
    if n <= 0:
        raise ValueError("Number of links n must be positive.")
    if l <= 0:
        raise ValueError("Link length l must be positive.")

    # The force is derived from the principles of statistical mechanics for an
    # isolated system (constant entropy).
    # 1. Total Entropy S = S_config + S_kinetic = constant
    #    S_config ~ -3*k_B*x^2 / (2*n*l^2)
    #    S_kinetic ~ n*k_B*ln(K)
    # 2. From dS=0, we find K(x) = E(0) * exp(3*x^2 / (2*n^2*l^2))
    # 3. Force f = -(dK/dx) at constant S
    #    f = -E(0) * d/dx[exp(3*x^2 / (2*n^2*l^2))]
    #    f = -E(0) * exp(3*x^2 / (2*n^2*l^2)) * (6x / (2*n^2*l^2))
    #    f = -E(0) * exp(3*x^2 / (2*n^2*l^2)) * (3x / (n^2*l^2))

    n_squared_l_squared = (n**2) * (l**2)
    exponent_term = (3 * x**2) / (2 * n_squared_l_squared)
    pre_factor = -(3 * x * E0) / n_squared_l_squared
    
    force = pre_factor * math.exp(exponent_term)
    
    return force

# Main part of the script
if __name__ == '__main__':
    # Print the symbolic formula for the user
    print("The derived force law f(x) for a thermally isolated polymer is:")
    print("f(x) = - (3 * x * E(0) / (n^2 * l^2)) * exp(3 * x^2 / (2 * n^2 * l^2))")
    print("-" * 50)
    print("where:")
    print("  f(x) = Force at extension x")
    print("  x    = End-to-end separation")
    print("  l    = Length of each segment")
    print("  n    = Number of segments")
    print("  E(0) = Kinetic energy at zero extension (x=0)")
    print("-" * 50)
    
    # --- Example Calculation ---
    # Define some example parameters for a hypothetical polymer
    # A large polymer with 1000 segments of length 1 unit.
    num_links = 1000
    link_length = 1.0  # arbitrary units of length
    # Initial kinetic energy (e.g., related to room temp, arbitrary units)
    energy_at_zero_extension = 414 # Proportional to k_B * T, ~ (300 K) in units of 1e-23 J
    # Calculate the force when stretched to 5% of its total contour length (n*l)
    extension = 0.05 * num_links * link_length

    print("Example calculation with the following parameters:")
    print(f"  Number of segments (n): {num_links}")
    print(f"  Segment length (l): {link_length}")
    print(f"  Extension (x): {extension}")
    print(f"  Initial Kinetic Energy (E(0)): {energy_at_zero_extension}")
    print("-" * 50)

    # Calculate and print the force
    calculated_force = calculate_isolated_polymer_force(
        x=extension,
        l=link_length,
        n=num_links,
        E0=energy_at_zero_extension
    )
    
    print(f"Calculated Restoring Force f({extension}): {calculated_force:.6f}")
    
    # Let's show how the force changes for a smaller extension for comparison
    small_extension = 1.0
    small_force = calculate_isolated_polymer_force(
        x=small_extension,
        l=link_length,
        n=num_links,
        E0=energy_at_zero_extension
    )
    print(f"Calculated Restoring Force f({small_extension}): {small_force:.6f}")
    # Small extension force using the linear approximation f ~ -3*x*E0/(n^2*l^2)
    approx_force = - (3 * small_extension * energy_at_zero_extension) / (num_links**2 * link_length**2)
    print(f"Linear Approximation at x=1.0:       {approx_force:.6f} (shows exponential is close to 1)")
    
