import math

def calculate_ground_state_energy():
    """
    Calculates the ground state energy for 4 hard-core bosons on a 7-site ring.
    """
    # System parameters
    N = 7
    num_photons = 4
    # We will calculate the energy in units of J by setting J=1
    J = 1.0

    # Step 1: Calculate all single-particle energy levels
    # E_m = -2*J*cos(2*pi*m/N)
    single_particle_energies = []
    for m in range(N):
        k_m = 2 * math.pi * m / N
        E_m = -2 * J * math.cos(k_m)
        single_particle_energies.append((E_m, m))

    # Step 2: Sort the energies to find the lowest ones
    single_particle_energies.sort()

    # Step 3: The ground state energy is the sum of the lowest 'num_photons' energies
    ground_state_energy = 0
    occupied_levels_m = []
    for i in range(num_photons):
        ground_state_energy += single_particle_energies[i][0]
        occupied_levels_m.append(single_particle_energies[i][1])
    
    occupied_levels_m.sort()

    # Step 4: Print the detailed calculation process
    print("In the limit U -> infinity, the system is equivalent to 4 non-interacting fermions on a 7-site ring.")
    print("The ground state energy is the sum of the 4 lowest single-particle energies.")
    print(f"The occupied single-particle states correspond to quantum numbers m = {occupied_levels_m}.")
    print("\nThe ground state energy E_gs is calculated as follows:")
    
    # Extract the individual cosine values for the equation
    c0 = math.cos(2 * math.pi * 0 / N)
    c1 = math.cos(2 * math.pi * 1 / N)
    c2 = math.cos(2 * math.pi * 2 / N)
    c6 = math.cos(2 * math.pi * 6 / N)
    
    print("\nE_gs = E_0 + E_1 + E_6 + E_2")
    print(f"E_gs = -2*J*cos(2*pi*0/7) - 2*J*cos(2*pi*1/7) - 2*J*cos(2*pi*6/7) - 2*J*cos(2*pi*2/7)")
    print(f"E_gs = -2*J * ({c0:.4f} + {c1:.4f} + {c6:.4f} + {c2:.4f})")
    
    factor = c0 + c1 + c6 + c2
    print("\nSince cos(2*pi*6/7) = cos(2*pi*1/7), this simplifies to:")
    print("E_gs = -2*J * (cos(0) + 2*cos(2*pi/7) + cos(4*pi/7))")
    print(f"E_gs = -2*J * ({c0:.4f} + 2*{c1:.4f} + {c2:.4f})")
    print(f"E_gs = -2*J * ({factor:.4f})")
    
    # Final answer in terms of J
    final_coefficient = ground_state_energy / J
    print(f"\nThe final ground state energy is {final_coefficient:.4f} * J.")

if __name__ == '__main__':
    calculate_ground_state_energy()
