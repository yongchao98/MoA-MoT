import collections

def find_missing_elements():
    """
    This function identifies which chemical elements, known in antiquity,
    are not mentioned in Homer's Odyssey.
    """
    
    # Step 1 & 2: Define the set of elements known in antiquity,
    # excluding zinc and antimony as requested.
    # These elements were generally recognized in their elemental form.
    elements_known_in_antiquity = {
        "Carbon",    # Known as charcoal/soot, but not seen as an element like metals were.
        "Copper",
        "Gold",
        "Silver",
        "Lead",
        "Tin",
        "Iron",
        "Sulfur",
        "Mercury"    # Unlikely to have been known to Homeric Greeks.
    }

    # Step 3: Define the set of elements mentioned in the Odyssey.
    # Textual analysis of the epic reveals mentions of the following:
    elements_in_odyssey = {
        "Gold",      # Greek: chrysos
        "Silver",    # Greek: argyros
        "Copper",    # Greek: chalkos (also used for its alloy, bronze)
        "Tin",       # Greek: kassiteros
        "Iron",      # Greek: sidēros
        "Lead",      # Greek: molybdos
        "Sulfur"     # Greek: theion (as brimstone for purification)
    }

    # Step 4: Calculate which elements are in the first set but not the second.
    missing_elements = sorted(list(elements_known_in_antiquity.difference(elements_in_odyssey)))

    # Step 5: Print the final answer.
    print("The anciently known elements (excluding zinc and antimony) not mentioned in the Odyssey are:")
    for element in missing_elements:
        print(f"- {element}")

find_missing_elements()
<<<['Carbon', 'Mercury']>>>