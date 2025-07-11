def find_unmentioned_elements():
    """
    Identifies which of the anciently known chemical elements are not mentioned
    in Homer's Odyssey.
    """

    # Step 1: Define the set of nine elements known in their elemental form in antiquity.
    # Zinc and antimony are excluded as per the prompt.
    ancient_elements = {
        "Gold", "Silver", "Copper", "Iron",
        "Lead", "Tin", "Mercury", "Sulfur", "Carbon"
    }

    # Step 2: Define the set of elements mentioned in the Odyssey.
    # This includes metals and sulfur (θεεῖον/theeion, or brimstone).
    # While Carbon (as charcoal/soot) was used, it wasn't typically described
    # as a distinct substance in the same vein as metals or sulfur.
    odyssey_elements = {
        "Gold", "Silver", "Copper", "Iron", "Lead", "Tin", "Sulfur"
    }

    # Step 3: Find the difference between the two sets.
    # These are the elements known in antiquity but not mentioned in the poem.
    unmentioned_elements = ancient_elements.difference(odyssey_elements)

    # Step 4: Print the results in a clear format.
    print("Elements known in antiquity (excluding zinc and antimony):")
    print(", ".join(sorted(list(ancient_elements))))
    print("\nElements mentioned as distinct substances in the Odyssey:")
    print(", ".join(sorted(list(odyssey_elements))))
    print("\nElements from the ancient list NOT mentioned in the Odyssey:")
    
    # The final 'equation' is the set difference, and we output each member.
    if not unmentioned_elements:
        print("None.")
    else:
        for element in sorted(list(unmentioned_elements)):
            print(f"- {element}")

if __name__ == "__main__":
    find_unmentioned_elements()