import numpy as np

def calculate_ground_state_energy():
    """
    Calculates the ground state energy for a system of hard-core bosons on a ring.
    """
    N = 7  # Number of cavities
    M = 4  # Number of photons

    print(f"Calculating the ground state energy for N={N} sites and M={M} photons in the U -> infinity limit.")
    print("This is equivalent to finding the ground state of M non-interacting fermions on a ring.")
    print("-" * 50)

    # The single-particle energy levels are E_m = -2J * cos(2*pi*m/N).
    # We list the unique energy levels and their degeneracies.
    # For N=7, the distinct m values can be chosen as m = 0, 1, 2, 3.
    # m=0 is non-degenerate. m=1,2,3 correspond to degenerate pairs (m and N-m).
    print("The single-particle energy levels are E_m = -2J * cos(2*pi*m/N):")
    
    energies = {}
    # We use m = 0, 1, -1, 2, -2, ... for clarity
    m_values_for_levels = [0, 1, 2, 3]
    for m in m_values_for_levels:
        energy_coeff = -2 * np.cos(2 * np.pi * m / N)
        energies[m] = energy_coeff
        degeneracy = "non-degenerate" if m == 0 else "doubly-degenerate"
        m_str = f"m={m}" if m==0 else f"m=+/_{m}"
        print(f"  Level E_{m} ({m_str}, {degeneracy}): {energy_coeff:.6f} * J")

    # To find the ground state for M particles, we fill the lowest energy levels.
    # We create a list of all N single-particle orbitals and sort them by energy.
    all_orbitals = []
    for m in range(N):
        energy_coeff = -2 * np.cos(2 * np.pi * m / N)
        # Store energy and the corresponding m value
        all_orbitals.append({'energy_coeff': energy_coeff, 'm': m})
    
    # Sort orbitals by energy in ascending order
    sorted_orbitals = sorted(all_orbitals, key=lambda x: x['energy_coeff'])

    # The M lowest energy orbitals for the ground state
    ground_state_orbitals = sorted_orbitals[:M]
    
    print("\n" + "-" * 50)
    print(f"To find the ground state for M={M} photons, we fill the {M} lowest energy orbitals:")
    for i, orbital in enumerate(ground_state_orbitals):
        print(f"  {i+1}. Occupy orbital m={orbital['m']}: Energy = {orbital['energy_coeff']:.6f} * J")

    # Sum the energies of these orbitals
    total_energy_J_coeff = sum(o['energy_coeff'] for o in ground_state_orbitals)

    print("\n" + "-" * 50)
    print("The total energy from the hopping term (E_hopping) is the sum of these energies.")
    
    occupied_energies_str = " + ".join([f"({o['energy_coeff']:.6f})" for o in ground_state_orbitals])
    print(f"E_hopping = ({occupied_energies_str}) * J")
    print(f"E_hopping = {total_energy_J_coeff:.6f} * J")

    print("\nThe full Hamiltonian includes the on-site energy term sum(omega * a_i.dag * a_i),")
    print(f"which adds a constant offset of M * omega = {M} * \u03C9.")
    
    print("\n" + "=" * 50)
    print("The final ground state energy E_ground is:")
    print(f"E_ground = {M} * \u03C9 + E_hopping")
    print(f"E_ground = {M} * \u03C9 + ({occupied_energies_str}) * J")
    print(f"E_ground = {M} * \u03C9 + ({total_energy_J_coeff:.6f}) * J")
    print("=" * 50)

calculate_ground_state_energy()