def find_optimum_material():
    """
    Compares materials based on their radiation length to find the optimum
    choice for a particle detector cooling system.
    """
    # Data for the materials. Radiation length (X0) is a measure of how
    # far a high-energy electron travels before losing all but 1/e of its energy.
    # A longer radiation length is better for a particle detector as it means
    # less interaction with the particles being measured.
    # Values are in centimeters (cm).
    materials = {
        'Titanium': 3.56,
        'Aluminium': 8.9,
        '316 Stainless Steel': 1.74, # Approximate value, similar to Iron
        'Copper': 1.43,
        'Nickel': 1.42
    }

    optimum_material_name = ''
    max_radiation_length = 0

    print("Comparing materials based on Radiation Length (X₀):")
    print("A longer radiation length is better to minimize particle interaction.\n")

    for name, length in materials.items():
        print(f"Material: {name}, Radiation Length: {length} cm")
        if length > max_radiation_length:
            max_radiation_length = length
            optimum_material_name = name

    print("\n--- Conclusion ---")
    print(f"The optimum choice is the material with the maximum radiation length.")
    print(f"Optimum Material: {optimum_material_name}")
    print(f"Maximum Radiation Length: {max_radiation_length} cm")

if __name__ == '__main__':
    find_optimum_material()