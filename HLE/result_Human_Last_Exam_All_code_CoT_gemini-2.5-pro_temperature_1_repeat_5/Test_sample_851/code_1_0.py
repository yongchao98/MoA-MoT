import pandas as pd

def find_optically_active_achiral_nonpolar_classes():
    """
    Identifies crystal point groups that are achiral, non-polar, and can exhibit optical activity.
    """
    # Data for the 32 crystallographic point groups
    # Properties:
    # 'chiral': Can exist in left- and right-handed forms. Lacks any improper rotation (m, -n).
    # 'polar': Has a unique vector direction (dipole moment is possible).
    # 'centrosymmetric': Possesses a center of inversion (-1).
    point_groups_data = {
        # Triclinic
        '1':    {'chiral': True,  'polar': True,  'centrosymmetric': False},
        '-1':   {'chiral': False, 'polar': False, 'centrosymmetric': True},
        # Monoclinic
        '2':    {'chiral': True,  'polar': True,  'centrosymmetric': False},
        'm':    {'chiral': False, 'polar': True,  'centrosymmetric': False},
        '2/m':  {'chiral': False, 'polar': False, 'centrosymmetric': True},
        # Orthorhombic
        '222':  {'chiral': True,  'polar': False, 'centrosymmetric': False},
        'mm2':  {'chiral': False, 'polar': True,  'centrosymmetric': False},
        'mmm':  {'chiral': False, 'polar': False, 'centrosymmetric': True},
        # Tetragonal
        '4':    {'chiral': True,  'polar': True,  'centrosymmetric': False},
        '-4':   {'chiral': False, 'polar': False, 'centrosymmetric': False},
        '4/m':  {'chiral': False, 'polar': False, 'centrosymmetric': True},
        '422':  {'chiral': True,  'polar': False, 'centrosymmetric': False},
        '4mm':  {'chiral': False, 'polar': True,  'centrosymmetric': False},
        '-42m': {'chiral': False, 'polar': False, 'centrosymmetric': False},
        '4/mmm':{'chiral': False, 'polar': False, 'centrosymmetric': True},
        # Trigonal
        '3':    {'chiral': True,  'polar': True,  'centrosymmetric': False},
        '-3':   {'chiral': False, 'polar': False, 'centrosymmetric': True},
        '32':   {'chiral': True,  'polar': False, 'centrosymmetric': False},
        '3m':   {'chiral': False, 'polar': True,  'centrosymmetric': False},
        '-3m':  {'chiral': False, 'polar': False, 'centrosymmetric': True},
        # Hexagonal
        '6':    {'chiral': True,  'polar': True,  'centrosymmetric': False},
        '-6':   {'chiral': False, 'polar': False, 'centrosymmetric': False},
        '6/m':  {'chiral': False, 'polar': False, 'centrosymmetric': True},
        '622':  {'chiral': True,  'polar': False, 'centrosymmetric': False},
        '6mm':  {'chiral': False, 'polar': True,  'centrosymmetric': False},
        '-6m2': {'chiral': False, 'polar': False, 'centrosymmetric': False}, # Also written as -62m
        '6/mmm':{'chiral': False, 'polar': False, 'centrosymmetric': True},
        # Cubic
        '23':   {'chiral': True,  'polar': False, 'centrosymmetric': False},
        'm-3':  {'chiral': False, 'polar': False, 'centrosymmetric': True},
        '432':  {'chiral': True,  'polar': False, 'centrosymmetric': False},
        '-43m': {'chiral': False, 'polar': False, 'centrosymmetric': False},
        'm-3m': {'chiral': False, 'polar': False, 'centrosymmetric': True},
    }

    print("Finding crystal classes based on the following criteria:")
    print("1. Achiral: The crystal is superimposable on its mirror image (property 'chiral' is False).")
    print("2. Non-polar: The crystal has no unique polar axis (property 'polar' is False).")
    print("3. Allows Optical Activity: The crystal is non-centrosymmetric (property 'centrosymmetric' is False).\n")

    qualified_classes = []
    for name, properties in point_groups_data.items():
        # Apply the three conditions
        is_achiral = not properties['chiral']
        is_non_polar = not properties['polar']
        is_optically_active = not properties['centrosymmetric']

        if is_achiral and is_non_polar and is_optically_active:
            qualified_classes.append(name)
    
    print("The crystal classes that are achiral, non-polar, and can be optically active are:")
    print(', '.join(sorted(qualified_classes)))
    print("\nNow, let's evaluate the given answer choices against this list.")
    
    options = {
        'A': ['m', 'mm2'],
        'B': ['-6', '-6m2', '-43m'],
        'C': ['3m', '4mm', '6mm'], # Note: 4m is commonly 4mm
        'D': ['-4', '-42m'],
        'E': ['1', '2', '3', '4', '6']
    }

    print("\nAnalysis of options:")
    for opt, classes in options.items():
        is_correct = all(c in qualified_classes for c in classes)
        # Check properties for explanation
        props = []
        for c in classes:
            p = point_groups_data.get(c, {})
            prop_str = f"Chiral: {p.get('chiral', 'N/A')}, Polar: {p.get('polar', 'N/A')}"
            props.append(prop_str)
        
        # Check against our derived list
        in_list_check = [c in qualified_classes for c in classes]

        print(f"Option {opt}: {classes}")
        if all(in_list_check):
             print(" -> This option is CORRECT. All classes are achiral, non-polar, and non-centrosymmetric.")
        else:
             print(" -> This option is INCORRECT.")
             for i, c in enumerate(classes):
                 if not in_list_check[i]:
                    p = point_groups_data.get(c, {})
                    if p.get('chiral'):
                        print(f"    - '{c}' is chiral.")
                    elif p.get('polar'):
                        print(f"    - '{c}' is polar.")
                    elif p.get('centrosymmetric'):
                         print(f"    - '{c}' is centrosymmetric and thus not optically active.")


if __name__ == '__main__':
    find_optically_active_achiral_nonpolar_classes()
