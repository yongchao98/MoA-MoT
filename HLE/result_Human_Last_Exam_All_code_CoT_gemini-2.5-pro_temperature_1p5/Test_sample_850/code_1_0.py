import collections

def find_crystal_classes():
    """
    Analyzes the 32 crystal point groups to find which are achiral, non-polar,
    and optically active.

    This function is based on the fundamental principles of crystallography:
    1.  Polarity: A crystal class is polar if it has a unique vector direction
        that is not cancelled by other symmetry operations. The polar classes
        are 1, 2, 3, 4, 6, m, mm2, 3m, 4mm, 6mm.
    2.  Chirality & Optical Activity: A crystal is chiral if its structure cannot
        be superimposed on its mirror image. This occurs in classes that lack
        any roto-inversion axes (e.g., a center of inversion 'i' or a mirror
        plane 'm'). Optical activity (the ability to rotate polarized light)
        is a direct consequence of chirality. Therefore, only chiral crystal
        classes can be optically active. An achiral crystal CANNOT be
        optically active.

    This script will demonstrate that the set of crystal classes requested is empty.
    """
    
    # Define the 32 crystal classes and their properties.
    # Key: Hermann-Mauguin symbol
    # Value: A dictionary with boolean properties for 'is_chiral' and 'is_polar'.
    # 'is_optically_active' is determined directly from 'is_chiral'.
    
    crystal_classes = {
        # Triclinic
        "1":    {'is_chiral': True,  'is_polar': True},
        "-1":   {'is_chiral': False, 'is_polar': False},
        # Monoclinic
        "2":    {'is_chiral': True,  'is_polar': True},
        "m":    {'is_chiral': False, 'is_polar': True},
        "2/m":  {'is_chiral': False, 'is_polar': False},
        # Orthorhombic
        "222":  {'is_chiral': True,  'is_polar': False},
        "mm2":  {'is_chiral': False, 'is_polar': True},
        "mmm":  {'is_chiral': False, 'is_polar': False},
        # Tetragonal
        "4":    {'is_chiral': True,  'is_polar': True},
        "-4":   {'is_chiral': False, 'is_polar': False},
        "4/m":  {'is_chiral': False, 'is_polar': False},
        "422":  {'is_chiral': True,  'is_polar': False},
        "4mm":  {'is_chiral': False, 'is_polar': True},
        "-42m": {'is_chiral': False, 'is_polar': False},
        "4/mmm":{'is_chiral': False, 'is_polar': False},
        # Trigonal
        "3":    {'is_chiral': True,  'is_polar': True},
        "-3":   {'is_chiral': False, 'is_polar': False},
        "32":   {'is_chiral': True,  'is_polar': False},
        "3m":   {'is_chiral': False, 'is_polar': True},
        "-3m":  {'is_chiral': False, 'is_polar': False},
        # Hexagonal
        "6":    {'is_chiral': True,  'is_polar': True},
        "-6":   {'is_chiral': False, 'is_polar': False},
        "6/m":  {'is_chiral': False, 'is_polar': False},
        "622":  {'is_chiral': True,  'is_polar': False},
        "6mm":  {'is_chiral': False, 'is_polar': True},
        "-6m2": {'is_chiral': False, 'is_polar': False},
        "6/mmm":{'is_chiral': False, 'is_polar': False},
        # Cubic
        "23":   {'is_chiral': True,  'is_polar': False},
        "m-3":  {'is_chiral': False, 'is_polar': False},
        "432":  {'is_chiral': True,  'is_polar': False},
        "-43m": {'is_chiral': False, 'is_polar': False},
        "m-3m": {'is_chiral': False, 'is_polar': False},
    }

    print("Analyzing the 32 crystal point groups based on symmetry rules.\n")
    print("Rule: A crystal class is optically active if and only if it is chiral.")
    print("This means any achiral class cannot be optically active.\n")
    print("Searching for crystal classes with the following properties:")
    print(" - Achiral (is_chiral = False)")
    print(" - Non-polar (is_polar = False)")
    print(" - Optically Active (is_optically_active = True)\n")

    found_classes = []
    
    for name, properties in crystal_classes.items():
        is_chiral = properties['is_chiral']
        is_polar = properties['is_polar']
        
        # The key rule: Optical activity requires chirality.
        is_optically_active = is_chiral
        
        # Check if the class matches the user's criteria.
        # This will never be true because if is_chiral is False, is_optically_active is also False.
        if not is_chiral and not is_polar and is_optically_active:
            found_classes.append(name)
            
    if not found_classes:
        print("Result: No crystal classes were found that are simultaneously achiral, non-polar, and optically active.")
        print("This is because the condition for optical activity is chirality. An achiral crystal cannot be optically active.")
    else:
        # This part of the code is unreachable due to the physical principles.
        print(f"Found matching crystal classes: {', '.join(found_classes)}")
        
if __name__ == "__main__":
    find_crystal_classes()
