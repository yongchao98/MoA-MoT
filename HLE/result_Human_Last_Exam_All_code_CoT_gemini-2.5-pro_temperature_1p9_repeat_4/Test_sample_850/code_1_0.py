import collections

def find_special_crystal_classes():
    """
    Identifies crystal classes that are achiral, non-polar, and can exhibit optical activity.

    The 32 crystallographic point groups are categorized based on their properties:
    - Optical Activity: Possible only in non-centrosymmetric classes (lacking a center of inversion, i).
    - Achiral: The point group contains an inversion center (i) or a mirror plane (m), or other rotoinversion axes.
    - Non-polar: The class is not one of the 10 polar classes {1, 2, m, mm2, 3, 3m, 4, 4mm, 6, 6mm}.

    This script defines all 32 classes with these properties and filters them
    to find the ones that meet all three conditions.
    """

    # Define the properties for each of the 32 crystal point groups.
    # The boolean properties are: is_centrosymmetric, is_chiral, is_polar.
    ClassProperties = collections.namedtuple('ClassProperties', ['is_centrosymmetric', 'is_chiral', 'is_polar'])
    
    crystal_classes = {
        # Triclinic
        '1': ClassProperties(is_centrosymmetric=False, is_chiral=True, is_polar=True),
        '-1': ClassProperties(is_centrosymmetric=True, is_chiral=False, is_polar=False),
        # Monoclinic
        '2': ClassProperties(is_centrosymmetric=False, is_chiral=True, is_polar=True),
        'm': ClassProperties(is_centrosymmetric=False, is_chiral=False, is_polar=True),
        '2/m': ClassProperties(is_centrosymmetric=True, is_chiral=False, is_polar=False),
        # Orthorhombic
        '222': ClassProperties(is_centrosymmetric=False, is_chiral=True, is_polar=False),
        'mm2': ClassProperties(is_centrosymmetric=False, is_chiral=False, is_polar=True),
        'mmm': ClassProperties(is_centrosymmetric=True, is_chiral=False, is_polar=False),
        # Tetragonal
        '4': ClassProperties(is_centrosymmetric=False, is_chiral=True, is_polar=True),
        '-4': ClassProperties(is_centrosymmetric=False, is_chiral=False, is_polar=False),
        '4/m': ClassProperties(is_centrosymmetric=True, is_chiral=False, is_polar=False),
        '422': ClassProperties(is_centrosymmetric=False, is_chiral=True, is_polar=False),
        '4mm': ClassProperties(is_centrosymmetric=False, is_chiral=False, is_polar=True),
        '-42m': ClassProperties(is_centrosymmetric=False, is_chiral=False, is_polar=False),
        '4/mmm': ClassProperties(is_centrosymmetric=True, is_chiral=False, is_polar=False),
        # Trigonal
        '3': ClassProperties(is_centrosymmetric=False, is_chiral=True, is_polar=True),
        '-3': ClassProperties(is_centrosymmetric=True, is_chiral=False, is_polar=False),
        '32': ClassProperties(is_centrosymmetric=False, is_chiral=True, is_polar=False),
        '3m': ClassProperties(is_centrosymmetric=False, is_chiral=False, is_polar=True),
        '-3m': ClassProperties(is_centrosymmetric=True, is_chiral=False, is_polar=False),
        # Hexagonal
        '6': ClassProperties(is_centrosymmetric=False, is_chiral=True, is_polar=True),
        '-6': ClassProperties(is_centrosymmetric=False, is_chiral=False, is_polar=False),
        '6/m': ClassProperties(is_centrosymmetric=True, is_chiral=False, is_polar=False),
        '622': ClassProperties(is_centrosymmetric=False, is_chiral=True, is_polar=False),
        '6mm': ClassProperties(is_centrosymmetric=False, is_chiral=False, is_polar=True),
        '-6m2': ClassProperties(is_centrosymmetric=False, is_chiral=False, is_polar=False),
        '6/mmm': ClassProperties(is_centrosymmetric=True, is_chiral=False, is_polar=False),
        # Cubic
        '23': ClassProperties(is_centrosymmetric=False, is_chiral=True, is_polar=False),
        'm-3': ClassProperties(is_centrosymmetric=True, is_chiral=False, is_polar=False),
        '432': ClassProperties(is_centrosymmetric=False, is_chiral=True, is_polar=False),
        '-43m': ClassProperties(is_centrosymmetric=False, is_chiral=False, is_polar=False),
        'm-3m': ClassProperties(is_centrosymmetric=True, is_chiral=False, is_polar=False),
    }

    result_classes = []
    for name, properties in crystal_classes.items():
        # Condition 1: Achiral (NOT chiral)
        is_achiral = not properties.is_chiral
        
        # Condition 2: Non-polar
        is_non_polar = not properties.is_polar
        
        # Condition 3: Allows optical activity (NOT centrosymmetric)
        allows_optical_activity = not properties.is_centrosymmetric

        if is_achiral and is_non_polar and allows_optical_activity:
            result_classes.append(name)
            
    print("The achiral, non-polar crystal classes with the correct symmetry for optical activity are:")
    # Print each class on its own line for clarity
    for class_name in result_classes:
        print(class_name)

if __name__ == '__main__':
    find_special_crystal_classes()