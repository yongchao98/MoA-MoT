def find_crystal_classes():
    """
    Analyzes the 32 crystal classes to find those that are optically active
    and non-polar, while explaining the key physics concepts involved.
    """

    # Data for the 32 crystal classes in Hermann-Mauguin notation.
    # 'chiral' is a prerequisite for optical activity.
    # 'polar' refers to the presence of a unique polar axis.
    crystal_classes = {
        # System: Triclinic
        '1':     {'chiral': True,  'polar': True, 'system': 'Triclinic'},
        '-1':    {'chiral': False, 'polar': False, 'system': 'Triclinic'},
        # System: Monoclinic
        '2':     {'chiral': True,  'polar': True, 'system': 'Monoclinic'},
        'm':     {'chiral': False, 'polar': True, 'system': 'Monoclinic'},
        '2/m':   {'chiral': False, 'polar': False, 'system': 'Monoclinic'},
        # System: Orthorhombic
        '222':   {'chiral': True,  'polar': False, 'system': 'Orthorhombic'},
        'mm2':   {'chiral': False, 'polar': True, 'system': 'Orthorhombic'},
        'mmm':   {'chiral': False, 'polar': False, 'system': 'Orthorhombic'},
        # System: Tetragonal
        '4':     {'chiral': True,  'polar': True, 'system': 'Tetragonal'},
        '-4':    {'chiral': False, 'polar': False, 'system': 'Tetragonal'},
        '4/m':   {'chiral': False, 'polar': False, 'system': 'Tetragonal'},
        '422':   {'chiral': True,  'polar': False, 'system': 'Tetragonal'},
        '4mm':   {'chiral': False, 'polar': True, 'system': 'Tetragonal'},
        '-42m':  {'chiral': False, 'polar': False, 'system': 'Tetragonal'},
        '4/mmm': {'chiral': False, 'polar': False, 'system': 'Tetragonal'},
        # System: Trigonal
        '3':     {'chiral': True,  'polar': True, 'system': 'Trigonal'},
        '-3':    {'chiral': False, 'polar': False, 'system': 'Trigonal'},
        '32':    {'chiral': True,  'polar': False, 'system': 'Trigonal'},
        '3m':    {'chiral': False, 'polar': True, 'system': 'Trigonal'},
        '-3m':   {'chiral': False, 'polar': False, 'system': 'Trigonal'},
        # System: Hexagonal
        '6':     {'chiral': True,  'polar': True, 'system': 'Hexagonal'},
        '-6':    {'chiral': False, 'polar': False, 'system': 'Hexagonal'},
        '6/m':   {'chiral': False, 'polar': False, 'system': 'Hexagonal'},
        '622':   {'chiral': True,  'polar': False, 'system': 'Hexagonal'},
        '6mm':   {'chiral': False, 'polar': True, 'system': 'Hexagonal'},
        '-6m2':  {'chiral': False, 'polar': False, 'system': 'Hexagonal'},
        '6/mmm': {'chiral': False, 'polar': False, 'system': 'Hexagonal'},
        # System: Cubic
        '23':    {'chiral': True,  'polar': False, 'system': 'Cubic'},
        'm-3':   {'chiral': False, 'polar': False, 'system': 'Cubic'},
        '432':   {'chiral': True,  'polar': False, 'system': 'Cubic'},
        '-43m':  {'chiral': False, 'polar': False, 'system': 'Cubic'},
        'm-3m':  {'chiral': False, 'polar': False, 'system': 'Cubic'},
    }

    # --- Explanation ---
    print("### Explanation of Crystal Symmetry and Optical Activity ###")
    print("\nYour question asks for crystal classes that are achiral, non-polar, and optically active.")
    print("This combination is impossible due to the fundamental requirements for optical activity.")
    print("\n1.  **Optical Activity requires Chirality:** A crystal can only be optically active if it is chiral. A chiral structure lacks symmetry elements like mirror planes (m) and centers of inversion (i or -1).")
    print("2.  **Achiral means Not Chiral:** An achiral crystal, by definition, possesses at least one of these symmetry elements (m or i).")
    print("\nTherefore, a crystal class cannot be both achiral and optically active. The conditions are mutually exclusive. There are **zero** crystal classes that satisfy your original request.")
    print("\n-------------------------------------------------------------")

    # --- Answering the Likely Intended Question ---
    print("\nIt is likely the question intended to ask for **chiral** (and thus optically active) and **non-polar** crystal classes.")
    print("The script below finds these specific classes.")

    intended_classes = []
    for name, properties in crystal_classes.items():
        # A class is optically active if it is chiral.
        is_optically_active = properties['chiral']
        is_non_polar = not properties['polar']

        if is_optically_active and is_non_polar:
            intended_classes.append(name)

    print("\n### Optically Active (Chiral) and Non-Polar Crystal Classes ###")
    print("The following crystal classes have the correct symmetry to be both optically active and non-polar:")
    print("---------------------------------------------------------------------")
    # Output each class name, interpreted as "each number in the final equation"
    for crystal_class in intended_classes:
        system = crystal_classes[crystal_class]['system']
        print(f"Class: {crystal_class:<5} (System: {system})")
    print("---------------------------------------------------------------------")

if __name__ == '__main__':
    find_crystal_classes()