def find_unmentioned_elements():
    """
    This function determines which elements known in antiquity are not mentioned in the Odyssey.
    """
    # Step 1: Define the set of elements known in their elemental form around the 8th-7th century BCE.
    # Zinc and Antimony are excluded as per the prompt.
    known_elements = {"Gold", "Silver", "Copper", "Iron", "Lead", "Tin", "Mercury", "Sulfur", "Carbon"}

    # Step 2: Define the set of elements mentioned in the Odyssey.
    # Gold (χρυσός), Silver (ἄργυρος), Copper (χαλκός), Iron (σίδηρος),
    # Tin (κασσίτερος), Sulfur (θεῖον), and Carbon (ἄνθραξ) are mentioned.
    mentioned_in_odyssey = {"Gold", "Silver", "Copper", "Iron", "Tin", "Sulfur", "Carbon"}

    # Step 3: Find the elements in the 'known' set but not in the 'mentioned' set.
    unmentioned_elements = known_elements.difference(mentioned_in_odyssey)

    # Step 4: Print the results clearly.
    print("Of the elements known in antiquity, the ones not mentioned in the Odyssey are:")
    # Sort the list for consistent and readable output
    for element in sorted(list(unmentioned_elements)):
        print(f"- {element}")

if __name__ == "__main__":
    find_unmentioned_elements()