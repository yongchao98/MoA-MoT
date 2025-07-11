def find_unmentioned_elements():
    """
    Identifies which chemical elements known in Homeric times are not mentioned in the Odyssey.
    """
    # Step 1: Define elements known in their elemental form around the 8th-7th century BCE.
    # This list includes the seven "Metals of Antiquity" plus Carbon and Sulfur.
    # Zinc and Antimony are excluded as per the user's request.
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

    # Step 2: Define elements mentioned in the Odyssey.
    # Copper and Tin are mentioned via "bronze" (an alloy of the two).
    # Carbon is mentioned as "charcoal" for fires.
    # Sulfur is used for fumigation.
    mentioned_in_odyssey = {
        "Gold",      # Mentioned frequently (e.g., golden cups).
        "Silver",    # Mentioned frequently (e.g., silver bowls).
        "Copper",    # Mentioned as a component of bronze, the primary metal for weapons and armor.
        "Tin",       # Also a component of bronze.
        "Iron",      # Mentioned explicitly (e.g., the axe-head trial, trade goods).
        "Carbon",    # Mentioned as charcoal for fires and smithing.
        "Sulfur"     # Mentioned by Odysseus for purifying his hall after killing the suitors.
    }

    # Step 3: Find the difference between the two sets.
    unmentioned_elements = sorted(list(known_elements - mentioned_in_odyssey))

    # Step 4: Print the result.
    print("Elements known in antiquity (excluding Zinc and Antimony) but not mentioned in the Odyssey:")
    for element in unmentioned_elements:
        print(f"- {element}")

find_unmentioned_elements()