def find_dinosaur_chemical():
    """
    This function identifies the polycyclic hydrocarbon named after a creature
    extinct for over 65 million years.
    """

    # The creature is a dinosaur.
    creature_part = "Dinosaur"

    # The "-ane" suffix is standard for saturated hydrocarbons.
    suffix_part = "ane"

    # Combine the parts to form the chemical's common name.
    chemical_name = creature_part + suffix_part

    print(f"The creature part of the name is: {creature_part}")
    print(f"The hydrocarbon suffix part is: {suffix_part}")
    print(f"The full name of the polycyclic hydrocarbon is: {chemical_name}")

# Run the function to display the answer.
find_dinosaur_chemical()