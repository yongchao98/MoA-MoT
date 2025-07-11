def find_unmentioned_elements():
    """
    This script identifies which chemical elements, known in antiquity,
    are not mentioned in Homer's Odyssey.

    The logic is as follows:
    1. A set of elements known in their elemental form before/around the
       time of the Odyssey's composition is created.
       (Excludes Zinc and Antimony as they were not fully recognized as distinct elements).
    2. A set of elements that are mentioned in the Odyssey is created based
       on textual analysis of the epic poem.
    3. The script finds the difference between these two sets to identify the
       elements that were known but not mentioned.
    4. The result is printed.
    """

    # Elements known to humans in elemental form in antiquity (c. 8th century BC)
    ancient_elements = {"Gold", "Silver", "Copper", "Iron", "Lead", "Tin", "Mercury", "Sulfur", "Carbon"}

    # Elements from the above list that are mentioned in the Odyssey
    # Gold, Silver, Copper (as bronze), Iron, Tin, and Lead are metals for tools, treasure, and weapons.
    # Sulfur (as brimstone) is used for purification.
    # Carbon (as charcoal) is the fuel for forges and fires.
    mentioned_in_odyssey = {"Gold", "Silver", "Copper", "Iron", "Lead", "Tin", "Sulfur", "Carbon"}

    # Find the elements in ancient_elements that are NOT in mentioned_in_odyssey
    unmentioned_elements = ancient_elements.difference(mentioned_in_odyssey)

    print("List of anciently known elements not mentioned in the Odyssey:")
    if unmentioned_elements:
        for element in unmentioned_elements:
            print(f"- {element}")
    else:
        print("All anciently known elements are mentioned in the Odyssey.")

find_unmentioned_elements()