def calculate_fluorine_atoms_in_perfluoronanocar():
    """
    Calculates the number of fluorine atoms in a hypothetical perfluoronanocar.
    """
    # The chemical formula for the original nanocar is C278H30.
    # We are interested in the number of hydrogen atoms.
    num_hydrogens_in_nanocar = 30

    # The prefix "perfluoro-" means every hydrogen atom (H) is replaced by a fluorine atom (F).
    # This is a one-to-one replacement.
    num_fluorine_replaces_hydrogen = 1

    # Calculate the total number of fluorine atoms in the perfluorinated version.
    total_fluorine_atoms = num_hydrogens_in_nanocar * num_fluorine_replaces_hydrogen

    # Print the explanation and the final equation.
    print("The base molecule, a nanocar, has 30 hydrogen atoms.")
    print("In a 'perfluoro-' compound, each hydrogen atom is replaced by one fluorine atom.")
    print("The calculation is as follows:")
    print(f"{num_hydrogens_in_nanocar} (hydrogens) * {num_fluorine_replaces_hydrogen} (fluorine per hydrogen) = {total_fluorine_atoms} (total fluorines)")
    print(f"\nTherefore, a hypothetical perfluoronanocar would contain {total_fluorine_atoms} fluorine atoms.")

calculate_fluorine_atoms_in_perfluoronanocar()