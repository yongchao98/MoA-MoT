import operator

def find_optimal_material():
    """
    Finds the optimal material for a particle detector cooling system
    based on maximizing the radiation length.
    """
    # Data: Material name and its approximate radiation length (X0) in cm.
    # A longer radiation length is better as it means less interaction with particles.
    materials = {
        'Titanium': 3.56,
        'Aluminium': 8.9,
        '316 Stainless Steel': 1.76,
        'Copper': 1.43,
        'Nickel': 1.42
    }

    print("Comparing materials based on their radiation length (X0).")
    print("A longer radiation length is desirable to minimize particle interaction.\n")

    for material, x0 in materials.items():
        print(f"{material}: Radiation Length = {x0} cm")

    # Find the material with the highest radiation length
    # The max function with a key returns the key from the dictionary that has the maximum value.
    optimal_material = max(materials, key=materials.get)
    max_x0 = materials[optimal_material]

    print(f"\nConclusion:")
    print(f"The material with the longest radiation length is {optimal_material} with a value of {max_x0} cm.")
    print(f"Therefore, for pipes of similar dimensions and considering only this unique parameter, {optimal_material} is the optimum choice.")

find_optimal_material()