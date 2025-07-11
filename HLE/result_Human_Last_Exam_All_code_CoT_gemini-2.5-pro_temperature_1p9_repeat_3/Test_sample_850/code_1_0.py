import textwrap

def analyze_crystal_classes():
    """
    Analyzes the 32 crystal classes to find those that are achiral, 
    non-polar, and optically active, and explains the results.
    """
    # Data for the 32 crystallographic point groups (crystal classes).
    # 'chiral' is a necessary condition for optical activity.
    point_groups = {
        # Triclinic
        '1':    {'chiral': True,  'polar': True},
        'ī':    {'chiral': False, 'polar': False},
        # Monoclinic
        '2':    {'chiral': True,  'polar': True},
        'm':    {'chiral': False, 'polar': True},
        '2/m':  {'chiral': False, 'polar': False},
        # Orthorhombic
        '222':  {'chiral': True,  'polar': False},
        'mm2':  {'chiral': False, 'polar': True},
        'mmm':  {'chiral': False, 'polar': False},
        # Tetragonal
        '4':    {'chiral': True,  'polar': True},
        '4̄':    {'chiral': False, 'polar': False},
        '4/m':  {'chiral': False, 'polar': False},
        '422':  {'chiral': True,  'polar': False},
        '4mm':  {'chiral': False, 'polar': True},
        '4̄2m':  {'chiral': False, 'polar': False},
        '4/mmm':{'chiral': False, 'polar': False},
        # Trigonal
        '3':    {'chiral': True,  'polar': True},
        '3̄':    {'chiral': False, 'polar': False},
        '32':   {'chiral': True,  'polar': False},
        '3m':   {'chiral': False, 'polar': True},
        '3̄m':   {'chiral': False, 'polar': False},
        # Hexagonal
        '6':    {'chiral': True,  'polar': True},
        '6̄':    {'chiral': False, 'polar': False},
        '6/m':  {'chiral': False, 'polar': False},
        '622':  {'chiral': True,  'polar': False},
        '6mm':  {'chiral': False, 'polar': True},
        '6̄m2':  {'chiral': False, 'polar': False},
        '6/mmm':{'chiral': False, 'polar': False},
        # Cubic
        '23':   {'chiral': True,  'polar': False},
        'm3̄':   {'chiral': False, 'polar': False},
        '432':  {'chiral': True,  'polar': False},
        '4̄3m':  {'chiral': False, 'polar': False},
        'm3̄m':  {'chiral': False, 'polar': False},
    }

    # --- Introduction and Explanation ---
    print("--- Analysis of Crystal Class Properties ---")
    
    explanation = (
        "Optical activity is a phenomenon that can only occur in CHIRAL crystal structures. "
        "A chiral class lacks an inversion center (ī) and mirror planes (m). "
        "An ACHIRAL class, by definition, possesses at least one of these symmetry elements, which forbids optical activity. "
        "Therefore, the request for a crystal class that is both 'achiral' and 'optically active' is a logical contradiction. "
        "This script will programmatically verify that no such class exists."
    )
    print(textwrap.fill(explanation, width=80))
    print("\n" + "="*80)

    # --- Step 1: Find optically active classes ---
    print("\nStep 1: Identify all classes with the correct symmetry for optical activity.")
    print("These are the 11 CHIRAL classes:")
    optically_active_classes = [name for name, props in point_groups.items() if props['chiral']]
    print("-> ", optically_active_classes)

    # --- Step 2: Find achiral, non-polar classes ---
    print("\nStep 2: Identify all ACHIRAL and NON-POLAR classes.")
    achiral_non_polar_classes = [name for name, props in point_groups.items() if not props['chiral'] and not props['polar']]
    print("-> ", achiral_non_polar_classes)

    # --- Step 3: Find the intersection ---
    print("\nStep 3: Find the intersection of the two sets above to find the answer.")
    print("Answering: What classes are in BOTH List 1 AND List 2?")
    
    # Calculate the intersection of the two lists
    result_set = set(optically_active_classes).intersection(set(achiral_non_polar_classes))
    
    if not result_set:
        print("\n--- CONCLUSION ---")
        print("The intersection is empty. There are ZERO crystal classes that satisfy all the conditions.")
    else:
        # This case is not possible given the physics
        print("\n--- CONCLUSION ---")
        print("The following crystal classes satisfy all conditions:", sorted(list(result_set)))

if __name__ == '__main__':
    analyze_crystal_classes()