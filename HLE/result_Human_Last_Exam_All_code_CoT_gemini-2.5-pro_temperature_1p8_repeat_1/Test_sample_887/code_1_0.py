import math

def analyze_reactivity():
    """
    This script illustrates how vibrational excitation of a C-H bond
    affects its reactivity in a reaction with atomic fluorine.
    """

    # --- Step 1: Define a simplified energy model (arbitrary units) ---
    # The activation energy barrier for the reaction
    reaction_barrier = 100

    # Ground state energy of the C-H bond
    ch_bond_ground_energy = 10

    # Ground state energy of the C-D bond (slightly lower due to heavier mass)
    cd_bond_ground_energy = 8

    # Energy added to the C-H bond by the IR laser
    laser_energy_input = 75

    print("--- Simplified Reaction Model ---")
    print(f"Energy of Reaction Barrier: {reaction_barrier} units")
    print(f"Ground State Energy of C-H bond: {ch_bond_ground_energy} units")
    print(f"Ground State Energy of C-D bond: {cd_bond_ground_energy} units")
    print("-" * 35)

    # --- Step 2: Calculate energy needed for reaction without laser ---
    # This is the energy gap that must be overcome for a reaction to occur.
    # A larger gap means lower reactivity.
    energy_needed_ch_ground = reaction_barrier - ch_bond_ground_energy
    energy_needed_cd_ground = reaction_barrier - cd_bond_ground_energy

    print("--- Reactivity in Ground State (No Laser) ---")
    print("Equation for energy gap for C-H bond:")
    print(f"{reaction_barrier} - {ch_bond_ground_energy} = {energy_needed_ch_ground} units")
    print("Equation for energy gap for C-D bond:")
    print(f"{reaction_barrier} - {cd_bond_ground_energy} = {energy_needed_cd_ground} units")
    print("-" * 35)

    # --- Step 3: Excite the C-H bond with the laser ---
    ch_bond_excited_energy = ch_bond_ground_energy + laser_energy_input
    print("--- C-H Bond Excitation by IR Laser ---")
    print("Energy is added ONLY to the C-H bond.")
    print("Equation for excited C-H bond energy:")
    print(f"{ch_bond_ground_energy} + {laser_energy_input} = {ch_bond_excited_energy} units")
    print("-" * 35)


    # --- Step 4: Calculate energy needed for reaction with excited C-H ---
    energy_needed_ch_excited = reaction_barrier - ch_bond_excited_energy

    print("--- Reactivity with Excited C-H Bond ---")
    print("The energy gap for the excited C-H bond is now much smaller.")
    print("Equation for energy gap for EXCITED C-H bond:")
    print(f"{reaction_barrier} - {ch_bond_excited_energy} = {energy_needed_ch_excited} units")
    print(f"Energy gap for unexcited C-D bond remains: {energy_needed_cd_ground} units")
    print("-" * 35)

    # --- Step 5: Conclusion ---
    print("\n--- Conclusion ---")
    print(f"Energy gap to reaction for excited C-H: {energy_needed_ch_excited} units")
    print(f"Energy gap to reaction for unexcited C-D: {energy_needed_cd_ground} units")
    print("\nBecause the energy gap for the excited C-H bond is significantly smaller,")
    print("its reactivity is massively increased, leading to much faster bond cleavage.")

analyze_reactivity()