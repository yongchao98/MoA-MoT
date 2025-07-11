import math

def calculate_ground_state_energy():
    """
    Calculates the ground state energy for a system of 4 hard-core bosons
    on a 7-site ring.
    """
    
    # Parameters
    N = 7  # Number of cavities
    M = 4  # Number of photons

    print("Step 1: System Definition")
    print(f"The system consists of {M} photons in {N} optical cavities arranged in a ring.")
    print("In the limit of infinite on-site interaction U, no two photons can occupy the same cavity.")
    print("This simplifies the system to {M} non-interacting fermions on a {N}-site ring.")
    print("-" * 50)
    
    print("Step 2: Energy Equation Structure")
    print("The total ground state energy E_ground is the sum of two parts:")
    print("1. A constant energy offset from the on-site energy term: M * omega")
    print("2. A hopping energy term E_hop from particles moving between cavities.")
    print(f"E_ground = {M}*omega + E_hop")
    print("-" * 50)
    
    print("Step 3: Calculating Hopping Energy E_hop")
    print("E_hop is found by filling the single-particle energy levels of a tight-binding ring.")
    print(f"The single-particle energy levels are e_m = -2*J*cos(2*pi*m/{N}).")
    print(f"For the ground state with {M} particles, we fill the lowest {M} available energy states:")
    print(" - 1 particle occupies the m=0 state.")
    print(" - 2 particles occupy the degenerate m=1 and m=6 states.")
    print(" - 1 particle occupies one of the degenerate m=2 and m=5 states.")
    print("\nThe energies of these occupied states (in units of J) are:")
    
    # Calculate single particle energies in units of J
    e0_J = -2 * math.cos(0)
    e1_J = -2 * math.cos(2 * math.pi / N)
    e2_J = -2 * math.cos(4 * math.pi / N)

    print(f"Energy of m=0 state: e_0 = {e0_J:.4f}*J")
    print(f"Energy of m=1,6 states: e_1 = {e1_J:.4f}*J")
    print(f"Energy of m=2,5 states: e_2 = {e2_J:.4f}*J")
    
    # The total hopping energy is the sum: E_hop = e_0 + 2*e_1 + e_2
    E_hop_J = e0_J + 2 * e1_J + e2_J
    
    print("\nThe total hopping energy is the sum of these individual energies:")
    print(f"E_hop = e_0 + 2*e_1 + e_2 = ({e0_J:.4f} + 2*({e1_J:.4f}) + {e2_J:.4f})*J")
    print(f"E_hop = ({E_hop_J:.4f})*J")
    print("-" * 50)
    
    print("Step 4: Final Ground State Energy Equation")
    print("Combining the offset and hopping terms gives the final equation:")
    print(f"E_ground = {M}*omega + ({E_hop_J:.4f})*J")

calculate_ground_state_energy()