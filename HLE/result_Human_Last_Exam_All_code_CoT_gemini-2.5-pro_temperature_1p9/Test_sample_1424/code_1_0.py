def main():
    """
    This script demonstrates the functional dependence of the bare Green's function (G_0)
    on the single-particle energy eigenvalue (ϵ_k) by showing the formula for
    both particle and hole propagation.
    """

    # --- Define physical parameters for the demonstration ---
    omega = 5.5       # Example frequency in units of energy
    epsilon_F = 5.0     # The Fermi energy
    # For the formula, eta (η) is an infinitesimal. We will represent it symbolically.
    eta_symbol = "η"

    print("The general formula for the bare Green's function is:")
    print("G_0(k, ω) = 1 / (ω - ϵ_k + iη * sgn(ϵ_k - ϵ_F))")
    print("-" * 50)

    # --- Case 1: Particle Propagation (Energy eigenvalue ABOVE the Fermi energy) ---
    # Here, we choose an energy epsilon_k that is greater than the Fermi energy.
    epsilon_k_particle = 6.0

    print("Case 1: Particle Propagation (since ϵ_k > ϵ_F)")
    print(f"Given values: ω = {omega}, ϵ_k = {epsilon_k_particle}, ϵ_F = {epsilon_F}")
    print("The specific formula for a particle is: G_0(k, ω) = 1 / (ω - ϵ_k + iη)")
    
    # We construct and print the final equation with the numbers substituted in.
    print("\nResulting equation:")
    print(f"G_0(k, ω) = 1 / ({omega} - {epsilon_k_particle} + i{eta_symbol})")
    
    print("\n" + "-" * 50)

    # --- Case 2: Hole Propagation (Energy eigenvalue BELOW the Fermi energy) ---
    # Here, we choose an energy epsilon_k that is less than the Fermi energy.
    epsilon_k_hole = 4.0
    
    print("Case 2: Hole Propagation (since ϵ_k < ϵ_F)")
    print(f"Given values: ω = {omega}, ϵ_k = {epsilon_k_hole}, ϵ_F = {epsilon_F}")
    print("The specific formula for a hole is: G_0(k, ω) = 1 / (ω - ϵ_k - iη)")

    # We construct and print the final equation with the numbers substituted in.
    print("\nResulting equation:")
    print(f"G_0(k, ω) = 1 / ({omega} - {epsilon_k_hole} - i{eta_symbol})")


if __name__ == "__main__":
    main()