import textwrap

def find_crystal_classes():
    """
    Analyzes the 32 crystal classes to find those that are simultaneously
    achiral, non-polar, and optically active. The script demonstrates
    that no such classes exist due to a logical contradiction in the requirements.
    """
    
    # (Hermann-Mauguin symbol, is_optically_active/chiral, is_polar)
    # Data for the 32 point groups. A class is optically active if and only if it is chiral.
    crystal_classes_data = {
        # Triclinic
        '1':        {'chiral': True, 'polar': True},
        'bar(1)':   {'chiral': False, 'polar': False},
        # Monoclinic
        '2':        {'chiral': True, 'polar': True},
        'm':        {'chiral': False, 'polar': True},
        '2/m':      {'chiral': False, 'polar': False},
        # Orthorhombic
        '222':      {'chiral': True, 'polar': False},
        'mm2':      {'chiral': False, 'polar': True},
        'mmm':      {'chiral': False, 'polar': False},
        # Tetragonal
        '4':        {'chiral': True, 'polar': True},
        'bar(4)':   {'chiral': False, 'polar': False},
        '4/m':      {'chiral': False, 'polar': False},
        '422':      {'chiral': True, 'polar': False},
        '4mm':      {'chiral': False, 'polar': True},
        'bar(4)2m': {'chiral': False, 'polar': False},
        '4/mmm':    {'chiral': False, 'polar': False},
        # Trigonal
        '3':        {'chiral': True, 'polar': True},
        'bar(3)':   {'chiral': False, 'polar': False},
        '32':       {'chiral': True, 'polar': False},
        '3m':       {'chiral': False, 'polar': True},
        'bar(3)m':  {'chiral': False, 'polar': False},
        # Hexagonal
        '6':        {'chiral': True, 'polar': True},
        'bar(6)':   {'chiral': False, 'polar': False},
        '6/m':      {'chiral': False, 'polar': False},
        '622':      {'chiral': True, 'polar': False},
        '6mm':      {'chiral': False, 'polar': True},
        'bar(6)m2': {'chiral': False, 'polar': False},
        '6/mmm':    {'chiral': False, 'polar': False},
        # Cubic
        '23':       {'chiral': True, 'polar': False},
        'm bar(3)': {'chiral': False, 'polar': False},
        '432':      {'chiral': True, 'polar': False},
        'bar(4)3m': {'chiral': False, 'polar': False},
        'm bar(3)m':{'chiral': False, 'polar': False}
    }

    explanation = """
    This program solves for crystal classes with three specific properties:
    1.  Optically Active: Can rotate polarized light. This requires a CHIRAL point group, meaning it has NO improper rotation axes (like mirror planes or inversion centers).
    2.  Achiral: Is superimposable on its mirror image. This means the point group MUST HAVE an improper rotation axis.
    3.  Non-polar: Lacks a unique directional axis.

    Based on definitions 1 and 2, there is a fundamental contradiction. A crystal class cannot simultaneously HAVE and NOT HAVE an improper rotation axis. Therefore, no crystal class can be both optically active and achiral.
    """
    print(textwrap.dedent(explanation).strip())
    
    print("\n--- Verifying with Data ---")

    matching_classes = []
    for name, properties in crystal_classes_data.items():
        is_optically_active = properties['chiral']
        is_achiral = not properties['chiral']
        is_non_polar = not properties['polar']
        
        if is_optically_active and is_achiral and is_non_polar:
            matching_classes.append(name)
            
    final_count = len(matching_classes)
    
    print("\nResult of search for classes that are (Optically Active AND Achiral AND Non-polar):")
    
    print("\nFinal Equation:")
    print(f"Number of matching crystal classes = {final_count}")

if __name__ == "__main__":
    find_crystal_classes()