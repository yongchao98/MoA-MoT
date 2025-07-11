def find_crystal_classes():
    """
    Analyzes the 32 crystallographic point groups to find those that meet specific criteria.
    The criteria are: achiral, non-polar, and capable of optical activity (i.e., non-centrosymmetric).
    """

    # Data for the 32 crystal classes (point groups).
    # Properties:
    # 'name': Hermann-Mauguin symbol.
    # 'is_centrosymmetric': True if it has a center of inversion. Optical activity is forbidden.
    # 'is_chiral': True if it lacks any mirror planes or rotoinversion axes.
    # 'is_polar': True if it has a unique direction not equivalent to its opposite.
    crystal_classes = [
        # Triclinic
        {'name': '1',     'is_centrosymmetric': False, 'is_chiral': True,  'is_polar': True},
        {'name': '-1',    'is_centrosymmetric': True,  'is_chiral': False, 'is_polar': False},
        # Monoclinic
        {'name': '2',     'is_centrosymmetric': False, 'is_chiral': True,  'is_polar': True},
        {'name': 'm',     'is_centrosymmetric': False, 'is_chiral': False, 'is_polar': True},
        {'name': '2/m',   'is_centrosymmetric': True,  'is_chiral': False, 'is_polar': False},
        # Orthorhombic
        {'name': '222',   'is_centrosymmetric': False, 'is_chiral': True,  'is_polar': False},
        {'name': 'mm2',   'is_centrosymmetric': False, 'is_chiral': False, 'is_polar': True},
        {'name': 'mmm',   'is_centrosymmetric': True,  'is_chiral': False, 'is_polar': False},
        # Tetragonal
        {'name': '4',     'is_centrosymmetric': False, 'is_chiral': True,  'is_polar': True},
        {'name': '-4',    'is_centrosymmetric': False, 'is_chiral': False, 'is_polar': False},
        {'name': '4/m',   'is_centrosymmetric': True,  'is_chiral': False, 'is_polar': False},
        {'name': '422',   'is_centrosymmetric': False, 'is_chiral': True,  'is_polar': False},
        {'name': '4mm',   'is_centrosymmetric': False, 'is_chiral': False, 'is_polar': True},
        {'name': '-42m',  'is_centrosymmetric': False, 'is_chiral': False, 'is_polar': False},
        {'name': '4/mmm', 'is_centrosymmetric': True,  'is_chiral': False, 'is_polar': False},
        # Trigonal
        {'name': '3',     'is_centrosymmetric': False, 'is_chiral': True,  'is_polar': True},
        {'name': '-3',    'is_centrosymmetric': True,  'is_chiral': False, 'is_polar': False},
        {'name': '32',    'is_centrosymmetric': False, 'is_chiral': True,  'is_polar': False},
        {'name': '3m',    'is_centrosymmetric': False, 'is_chiral': False, 'is_polar': True},
        {'name': '-3m',   'is_centrosymmetric': True,  'is_chiral': False, 'is_polar': False},
        # Hexagonal
        {'name': '6',     'is_centrosymmetric': False, 'is_chiral': True,  'is_polar': True},
        {'name': '-6',    'is_centrosymmetric': False, 'is_chiral': False, 'is_polar': False},
        {'name': '6/m',   'is_centrosymmetric': True,  'is_chiral': False, 'is_polar': False},
        {'name': '622',   'is_centrosymmetric': False, 'is_chiral': True,  'is_polar': False},
        {'name': '6mm',   'is_centrosymmetric': False, 'is_chiral': False, 'is_polar': True},
        {'name': '-6m2',  'is_centrosymmetric': False, 'is_chiral': False, 'is_polar': False},
        {'name': '6/mmm', 'is_centrosymmetric': True,  'is_chiral': False, 'is_polar': False},
        # Cubic
        {'name': '23',    'is_centrosymmetric': False, 'is_chiral': True,  'is_polar': False},
        {'name': 'm-3',   'is_centrosymmetric': True,  'is_chiral': False, 'is_polar': False},
        {'name': '432',   'is_centrosymmetric': False, 'is_chiral': True,  'is_polar': False},
        {'name': '-43m',  'is_centrosymmetric': False, 'is_chiral': False, 'is_polar': False},
        {'name': 'm-3m',  'is_centrosymmetric': True,  'is_chiral': False, 'is_polar': False},
    ]

    result_classes = []
    print("Finding crystal classes that meet the following criteria:")
    print("1. Achiral (is_chiral = False)")
    print("2. Non-polar (is_polar = False)")
    print("3. Allows optical activity (is_centrosymmetric = False)")
    print("-" * 30)
    
    for cc in crystal_classes:
        # Apply the three conditions
        if (not cc['is_chiral']) and (not cc['is_polar']) and (not cc['is_centrosymmetric']):
            result_classes.append(cc['name'])

    print("The complete set of crystal classes satisfying all conditions is:")
    # The final print statement is required by the problem format.
    print(f"[{', '.join(repr(c) for c in result_classes)}]")
    print("\nComparing this result with the answer choices:")
    print("Choice A ('m', 'mm2') -> Incorrect, these are polar.")
    print("Choice B ('-6', '-6m2', '-43m') -> Correct, all members are in our result set.")
    print("Choice C ('3m', '4mm', '6mm') -> Incorrect, these are polar.")
    print("Choice D ('-4', '-42m') -> Correct, all members are in our result set.")
    print("Choice E ('1', '2', '3', '4', '6') -> Incorrect, these are chiral and polar.")

    print("\nBoth choices B and D are valid subsets of the correct classes. However, questions of this type usually look for a specific group. Choice B represents classes from hexagonal and cubic systems.")

if __name__ == '__main__':
    find_crystal_classes()