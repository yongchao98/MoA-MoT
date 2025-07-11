def calculate_fluorine_atoms():
    """
    Calculates the number of fluorine atoms in a hypothetical perfluoronanocar.

    The calculation is based on the structure of a real nanocar from the Tour group.
    A "perfluorinated" molecule has all its hydrogen atoms replaced by fluorine atoms.
    We need to count the total number of hydrogen atoms on the nanocar's chassis.
    """

    # 1. Hydrogens on the central biphenyl core of the chassis.
    # The core is 1,1'-biphenyl, substituted at positions 2, 2', 5, and 5'.
    # Hydrogen atoms remain at positions 3, 3', 4, 4', 6, and 6'.
    h_on_core = 6
    print(f"Number of hydrogens on the central biphenyl core = {h_on_core}")

    # 2. Hydrogens on the four phenylene groups that form part of the axles.
    # Each phenylene ring is substituted at two positions, leaving 4 hydrogens.
    num_axle_phenylene_groups = 4
    h_per_axle_phenylene = 4
    h_on_axle_phenylenes = num_axle_phenylene_groups * h_per_axle_phenylene
    print(f"Number of hydrogens on the four axle phenylene groups = {h_on_axle_phenylenes}")

    # 3. Hydrogens on the four tert-butyl groups, -C(CH3)3.
    # Each tert-butyl group has 3 methyl groups, and each methyl has 3 hydrogens (3 * 3 = 9).
    num_tert_butyl_groups = 4
    h_per_tert_butyl = 9
    h_on_tert_butyls = num_tert_butyl_groups * h_per_tert_butyl
    print(f"Number of hydrogens on the four tert-butyl groups = {h_on_tert_butyls}")

    # 4. Total hydrogens, which equals the number of fluorines in the perfluorinated version.
    total_fluorines = h_on_core + h_on_axle_phenylenes + h_on_tert_butyls

    # 5. Print the final equation and the result.
    print(f"\nThe total number of fluorine atoms in the perfluoronanocar is the sum of these counts.")
    print(f"Total Fluorine Atoms = {h_on_core} + {h_on_axle_phenylenes} + {h_on_tert_butyls} = {total_fluorines}")

if __name__ == "__main__":
    calculate_fluorine_atoms()