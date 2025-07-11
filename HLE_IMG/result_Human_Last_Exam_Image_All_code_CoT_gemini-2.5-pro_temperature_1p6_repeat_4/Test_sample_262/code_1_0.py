def find_unstable_complexes():
    """
    Identifies iridium complexes with shorter expected lifetimes based on ligand stability.

    The lifetime of these emitters is linked to their chemical stability.
    Fluorination of the cyclometalating phenylpyridine (C^N) ligands enhances stability.
    The unsubstituted ligand, 2-phenylpyridine ('ppy'), is known to be less stable
    and leads to shorter lifetimes in devices.

    This script checks which complexes contain this less stable ligand.
    """

    # Defining the cyclometalating (C^N) ligands for each complex.
    # The 'ppy' ligand is identified by its chemical name '2-phenylpyridine'.
    complex_ligands = {
        1: ['2-(2,4-difluorophenyl)pyridine', '2-phenylpyridine'],
        2: ['2-(2,4-difluorophenyl)pyridine', '2-(2,4-difluorophenyl)pyridine'],
        3: ['2-(2-fluorophenyl)pyridine', '2-phenylpyridine'],
        4: ['2-(2,4-difluorophenyl)pyridine', '2-(2-fluorophenyl)pyridine']
    }

    unstable_ligand_name = '2-phenylpyridine'
    short_lifetime_complexes = []

    for complex_number, ligands in complex_ligands.items():
        if unstable_ligand_name in ligands:
            short_lifetime_complexes.append(complex_number)

    print("The complexes expected to show shorter lifetimes are those containing the least stable '2-phenylpyridine' ligand.")
    print(f"Complex {short_lifetime_complexes[0]} contains this ligand.")
    print(f"Complex {short_lifetime_complexes[1]} contains this ligand.")
    print(f"Therefore, the final answer is [{short_lifetime_complexes[0]}, {short_lifetime_complexes[1]}].")

find_unstable_complexes()