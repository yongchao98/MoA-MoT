def simulate_chd3_reaction():
    """
    This script models the effect of C-H bond excitation on the reaction
    of deuterated methane (CHD3) with a fluorine atom.
    """
    # Define arbitrary energy units for our conceptual model
    reaction_barrier = 10.0
    ground_state_bond_energy = 3.0
    laser_excitation_energy = 8.0

    # --- Case 1: Ground State Molecule (no laser) ---
    print("--- Case 1: Ground State CHD3 (No Laser Excitation) ---")
    ch_bond_energy_ground = ground_state_bond_energy
    cd_bond_energy_ground = ground_state_bond_energy

    print(f"Energy in C-H bond: {ch_bond_energy_ground}")
    print(f"Energy in C-D bonds: {cd_bond_energy_ground}")
    print(f"Reaction Barrier: {reaction_barrier}")

    print("\nResult:")
    if ch_bond_energy_ground < reaction_barrier:
        print("The C-H bond does not have enough energy to react.")
    if cd_bond_energy_ground < reaction_barrier:
        print("The C-D bonds do not have enough energy to react.")
    print("-" * 50)

    # --- Case 2: C-H Bond is Vibrationally Excited by a Laser ---
    print("\n--- Case 2: Excited State CHD3 (C-H Bond Excited by Laser) ---")
    ch_bond_energy_excited = ground_state_bond_energy + laser_excitation_energy
    cd_bond_energy_excited = ground_state_bond_energy # C-D bonds are not excited

    print(f"Energy in C-H bond: {ground_state_bond_energy} + {laser_excitation_energy} = {ch_bond_energy_excited}")
    print(f"Energy in C-D bonds: {cd_bond_energy_excited}")
    print(f"Reaction Barrier: {reaction_barrier}")

    print("\nResult:")
    if ch_bond_energy_excited >= reaction_barrier:
        print("The excited C-H bond has enough energy to react and break!")
        print("This greatly enhances the likelihood of H atom removal.")
    else:
        print("The C-H bond still does not have enough energy to react.")

    if cd_bond_energy_excited >= reaction_barrier:
        print("The C-D bonds have enough energy to react.")
    else:
        print("The unexcited C-D bonds do not have enough energy to react.")
    print("-" * 50)

    print("\nConclusion:")
    print("Vibrational excitation of the C-H bond localizes energy in that bond,")
    print("making it much more reactive than the C-D bonds. This accelerates the")
    print("reaction by specifically enhancing the probability of H atom removal over D atoms.")

# Run the simulation
simulate_chd3_reaction()