import operator

def find_optimal_material():
    """
    Identifies the optimal material for a particle detector's cooling system
    by maximizing the radiation length.
    """
    # Data for radiation length (X0) for various materials.
    # A longer radiation length is better as it means less interaction with particles.
    materials = {
        'A. Titanium': 3.56,  # cm
        'B. Aluminium': 8.9,   # cm
        'C. 316 Stainless Steel': 1.76, # cm (approximated as Iron)
        'D. Copper': 1.43,   # cm
        'E. Nickle': 1.42,   # cm
    }

    print("Comparing materials based on their radiation length (longer is better):")
    for material, length in materials.items():
        print(f"- {material}: {length} cm")

    # Find the material with the highest radiation length
    optimal_choice = max(materials.items(), key=operator.itemgetter(1))
    
    print(f"\nThe material with the longest radiation length is '{optimal_choice[0]}' with a value of {optimal_choice[1]} cm.")
    print("Therefore, it is the optimum choice to minimize interference with particle detection.")

find_optimal_material()