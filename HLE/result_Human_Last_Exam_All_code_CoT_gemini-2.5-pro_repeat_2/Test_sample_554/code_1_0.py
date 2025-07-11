def calculate_fluorine_atoms_in_perfluoronanocar():
    """
    Calculates the number of fluorine atoms in a hypothetical perfluoronanocar.
    """
    # Step 1: Define the number of hydrogen atoms in a standard nanocar.
    # The nanocar's chassis and axle structure has the chemical formula C44H28O4.
    # The four C60 fullerene wheels have no hydrogen atoms.
    # Therefore, the entire molecule has 28 hydrogen atoms.
    num_hydrogen_atoms = 28

    # Step 2: Apply the "perfluoro-" rule.
    # The prefix "perfluoro-" means every hydrogen atom is replaced by a fluorine atom.
    # Thus, the number of fluorine atoms is equal to the number of hydrogen atoms.
    num_fluorine_atoms = num_hydrogen_atoms

    # Step 3: Display the reasoning and the final equation.
    print("A standard nanocar molecule contains a chassis and four fullerene wheels.")
    print(f"The chassis and axles part contains all the hydrogen atoms.")
    print(f"Number of hydrogen atoms in a nanocar = {num_hydrogen_atoms}")
    print("\nThe prefix 'perfluoro-' means all hydrogen atoms are replaced by fluorine atoms.")
    print(f"Therefore, the number of fluorine atoms in a perfluoronanocar is also {num_fluorine_atoms}.")
    print("\nFinal Equation:")
    print(f"Number of Fluorine Atoms = {num_fluorine_atoms}")

calculate_fluorine_atoms_in_perfluoronanocar()