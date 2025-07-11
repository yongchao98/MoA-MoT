def find_unmentioned_elements():
    """
    This script identifies the ancient chemical elements not mentioned in Homer's Odyssey.
    """

    # Step 1: Define the set of elements known in elemental form around the 8th-7th century BC.
    # Zinc and Antimony are excluded as per the problem description.
    # The list includes the seven metals of antiquity plus carbon and sulfur.
    known_elements = {
        "Gold",
        "Silver",
        "Copper",
        "Tin",
        "Lead",
        "Iron",
        "Mercury",
        "Carbon",
        "Sulfur"
    }

    # Step 2: Define the set of elements mentioned in the Odyssey, based on literary analysis.
    # Copper and Tin are mentioned, often as their alloy, bronze.
    # Carbon is mentioned as charcoal, and Sulfur is used for purification.
    # Lead is mentioned, though very rarely.
    mentioned_in_odyssey = {
        "Gold",
        "Silver",
        "Copper",
        "Tin",
        "Lead",
        "Iron",
        "Carbon",
        "Sulfur"
    }

    # Step 3: Calculate the set difference to find elements that were known but not mentioned.
    unmentioned_elements = known_elements.difference(mentioned_in_odyssey)

    # Step 4: Print the result.
    print("Of the elements known in antiquity (excluding Zinc and Antimony), the following are not mentioned in the Odyssey:")
    if not unmentioned_elements:
        print("None. All known elements from the period are mentioned.")
    else:
        for element in unmentioned_elements:
            print(f"- {element}")

find_unmentioned_elements()