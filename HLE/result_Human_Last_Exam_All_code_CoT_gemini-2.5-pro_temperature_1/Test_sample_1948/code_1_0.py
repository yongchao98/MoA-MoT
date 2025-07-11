def find_unmentioned_elements():
    """
    This script identifies the chemical elements known in antiquity (at the time of
    the Odyssey's composition) but not mentioned in the epic poem.
    """

    # Elements known and isolated in their elemental form around the 8th century BCE.
    # Zinc and Antimony are excluded as per the user's request.
    known_elements = {"Gold", "Silver", "Copper", "Iron", "Lead", "Tin", "Mercury", "Carbon", "Sulfur"}

    # Elements mentioned in The Odyssey, either directly or by clear implication (e.g., bronze implies copper and tin).
    mentioned_elements = {
        "Gold",      # Mentioned frequently as a sign of wealth (chrysos).
        "Silver",    # Mentioned as a precious metal for bowls and treasure (argyros).
        "Copper",    # A primary component of bronze (chalkos), used for weapons, armor, and tripods.
        "Tin",       # A necessary component of bronze (kassiteros), its use is implied with bronze.
        "Iron",      # Used for tools and weapons, central to the contest of the axes (sideros).
        "Carbon",    # Mentioned in the form of charcoal (anthrax) for forges and fires.
        "Sulfur"     # Mentioned as "brimstone" (theion) used by Odysseus to purify his hall.
    }

    # Find the elements that are in the 'known_elements' set but not in the 'mentioned_elements' set.
    unmentioned_elements = sorted(list(known_elements.difference(mentioned_elements)))

    print("Of the elements known in their pure form during the time of The Odyssey's composition (excluding zinc and antimony), the ones not mentioned in the poem are:")
    for element in unmentioned_elements:
        print(f"- {element}")

find_unmentioned_elements()
<<<Lead and Mercury>>>