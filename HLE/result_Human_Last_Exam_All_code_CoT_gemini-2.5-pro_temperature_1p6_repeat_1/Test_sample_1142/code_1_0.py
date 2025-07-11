import math

def get_electron_configuration(atomic_number):
    """Calculates the electron configuration based on the standard aufbau principle."""
    electrons = atomic_number
    config = {}
    # Use standard 3D aufbau order and capacities as per the prompt's rules
    aufbau_order = ['1s', '2s', '2p', '3s', '3p', '4s', '3d', '4p', '5s', '4d', '5p', '6s', '4f', '5d', '6p']
    subshell_capacities = {'s': 2, 'p': 6, 'd': 10, 'f': 14}
    
    config_str_list = []

    for subshell in aufbau_order:
        if electrons <= 0:
            break
        
        subshell_type = subshell[-1]
        capacity = subshell_capacities[subshell_type]
        
        num_electrons_in_subshell = min(electrons, capacity)
        config[subshell] = num_electrons_in_subshell
        config_str_list.append(f"{subshell}{num_electrons_in_subshell}")
        electrons -= num_electrons_in_subshell
        
    return config, " ".join(config_str_list)

def get_bonds_to_complete_subshell(config):
    """Calculates bonds needed to complete the highest-energy, partially-filled subshell."""
    subshell_capacities = {'s': 2, 'p': 6, 'd': 10, 'f': 14}
    
    if not config:
        return 0, None, 0

    # Per the aufbau principle, the highest-energy subshell is the last one being filled.
    highest_energy_subshell = list(config.keys())[-1]
    
    electrons_in_subshell = config[highest_energy_subshell]
    subshell_type = highest_energy_subshell[-1]
    capacity = subshell_capacities[subshell_type]
    
    # If the subshell is full, the atom is stable and forms 0 bonds.
    if electrons_in_subshell == capacity:
        return 0, highest_energy_subshell, electrons_in_subshell

    bonds_needed = capacity - electrons_in_subshell
    return bonds_needed, highest_energy_subshell, electrons_in_subshell

# --- Main analysis ---

# Define atoms by atomic number
Z_C = 6
Z_Ni = 28

# --- Step 1: Analyze Carbon's bonding preference ---
c_config, c_config_str = get_electron_configuration(Z_C)
c_bonds_needed, c_subshell, c_electrons = get_bonds_to_complete_subshell(c_config)
subshell_capacities = {'s': 2, 'p': 6, 'd': 10, 'f': 14}
c_subshell_capacity = subshell_capacities[c_subshell[-1]]


print(f"Analysis for Carbon (Z={Z_C}):")
print(f"Electron configuration is {c_config_str}.")
print(f"The highest-energy subshell ({c_subshell}) has {c_electrons} of a possible {c_subshell_capacity} electrons.")
print(f"To achieve a completed subshell, Carbon requires {c_bonds_needed} additional electrons.")
print(f"Therefore, Carbon atoms prefer to form {c_bonds_needed} covalent bonds.\n")

# --- Step 1: Analyze Nickel's bonding preference ---
ni_config, ni_config_str = get_electron_configuration(Z_Ni)
ni_bonds_needed, ni_subshell, ni_electrons = get_bonds_to_complete_subshell(ni_config)
ni_subshell_capacity = subshell_capacities[ni_subshell[-1]]

print(f"Analysis for Nickel (Z={Z_Ni}):")
print(f"Electron configuration is {ni_config_str}.")
print(f"The highest-energy subshell ({ni_subshell}) has {ni_electrons} of a possible {ni_subshell_capacity} electrons.")
print(f"To achieve a completed subshell, Nickel requires {ni_bonds_needed} additional electrons.")
print(f"Therefore, Nickel atoms prefer to form {ni_bonds_needed} covalent bonds.\n")

# --- Step 2 & 3: Resolve conflict and select structure ---
print("Crystal Structure Analysis for NiC:")
print(f"In a 1:1 crystal, C atoms wanting {c_bonds_needed} bonds and Ni atoms wanting {ni_bonds_needed} bonds creates a conflict.")
print("A stable, regular lattice can be formed if the atoms adopt a compromise coordination number, which is the average of their ideal values.")

compromise_bonds = (c_bonds_needed + ni_bonds_needed) / 2
print(f"Average bonds = ({c_bonds_needed} + {ni_bonds_needed}) / 2 = {int(compromise_bonds)}")
print(f"We look for a structure from the options where atoms have a degree of {int(compromise_bonds)}.\n")

print("Reviewing options with degree=3:")
print("C. tiling by octagons and squares (3)")
print("D. tiling by hexagons (3)")
print("A hexagonal lattice (graphene-like) has uniform 120-degree bond angles, which are ideal for the sp2-like hybridization needed to form 3 bonds.")
print("This is generally more energetically favorable than the mixed angles of an octagon-square tiling.")
print("Therefore, 'D. tiling by hexagons' is the most probable structure.\n")

# --- Step 4: Analyze Shear Strength ---
print("Shear Strength Isotropy Analysis:")
print("The predicted structure is a 2D hexagonal lattice.")
print("In such a lattice, the arrangement of atoms and the forces between them are not the same in every direction. Properties depend on the direction of applied stress.")
print("This means the material is anisotropic.")
print("Therefore, the crystal shear strength is NOT nearly isotropic.\n")

# --- Final Answer ---
chosen_option = "D"
isotropy_answer = "no"

print(f"<<<{chosen_option} {isotropy_answer}>>>")