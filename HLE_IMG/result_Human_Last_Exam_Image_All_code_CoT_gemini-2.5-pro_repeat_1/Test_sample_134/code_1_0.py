def name_molecule_and_calculate_formula():
    """
    This script identifies the molecule, names it, and calculates its chemical formula.
    """
    
    # 1. Define atomic composition of each building block.
    # para-phenylene group (-C6H4-)
    c_phenylene = 6
    h_phenylene = 4
    # vinylene group (-CH=CH-)
    c_vinylene = 2
    h_vinylene = 2
    # ethynylene group (-C≡C-)
    c_ethynylene = 2
    h_ethynylene = 0

    # 2. Define the composition of the main repeating unit.
    # The repeating unit is [para-phenylene-vinylene-para-phenylene-ethynylene].
    num_phenylene_per_unit = 2
    num_vinylene_per_unit = 1
    num_ethynylene_per_unit = 1

    # 3. Define the number of repeating units in the macrocycle.
    # By counting, we find there are 5 such units forming the ring.
    num_units_in_cycle = 5

    # 4. Calculate the chemical formula for one repeating unit.
    c_per_unit = (num_phenylene_per_unit * c_phenylene) + \
                 (num_vinylene_per_unit * c_vinylene) + \
                 (num_ethynylene_per_unit * c_ethynylene)
                 
    h_per_unit = (num_phenylene_per_unit * h_phenylene) + \
                 (num_vinylene_per_unit * h_vinylene) + \
                 (num_ethynylene_per_unit * h_ethynylene)

    # 5. Calculate the total chemical formula for the entire molecule.
    total_c = num_units_in_cycle * c_per_unit
    total_h = num_units_in_cycle * h_per_unit

    # 6. Print the name and the calculation details.
    molecule_name = "cyclo[5](para-phenylenevinylene-para-phenyleneethynylene)"
    
    print(f"Molecule Name: {molecule_name}\n")
    print("Chemical Formula Calculation:")
    print("The molecule is a cyclic pentamer (n=5) with a repeating unit of [p-phenylene-vinylene-p-phenylene-ethynylene].")
    print("-" * 40)
    print("Atoms per repeating unit:")
    print(f"  Carbon (C): ({num_phenylene_per_unit} * {c_phenylene}) + ({num_vinylene_per_unit} * {c_vinylene}) + ({num_ethynylene_per_unit} * {c_ethynylene}) = {c_per_unit}")
    print(f"  Hydrogen (H): ({num_phenylene_per_unit} * {h_phenylene}) + ({num_vinylene_per_unit} * {h_vinylene}) + ({num_ethynylene_per_unit} * {h_ethynylene}) = {h_per_unit}")
    print(f"  Formula of one unit: C{c_per_unit}H{h_per_unit}")
    print("-" * 40)
    print("Total atoms in the molecule (5 units):")
    print(f"  Total Carbon (C): {num_units_in_cycle} * {c_per_unit} = {total_c}")
    print(f"  Total Hydrogen (H): {num_units_in_cycle} * {h_per_unit} = {total_h}")
    print("-" * 40)
    print(f"Final Chemical Formula: C{total_c}H{total_h}")

if __name__ == "__main__":
    name_molecule_and_calculate_formula()