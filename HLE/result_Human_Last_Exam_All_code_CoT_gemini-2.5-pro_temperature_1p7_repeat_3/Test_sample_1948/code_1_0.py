import collections

def find_unmentioned_elements():
    """
    This script identifies which chemical elements, known in antiquity,
    were not mentioned in Homer's Odyssey.
    """

    # A set of the nine elements known in elemental form during the Greek Dark Ages.
    # Zinc and Antimony are excluded as per the prompt's instructions.
    elements_known_in_antiquity = {"Carbon", "Copper", "Gold", "Silver", "Lead", "Tin", "Iron", "Sulfur", "Mercury"}

    # A set of the elements found to be mentioned in the Odyssey.
    # This includes mentions of bronze (Copper/Tin alloy), brimstone (Sulfur),
    # and charcoal (Carbon).
    elements_mentioned_in_odyssey = {"Carbon", "Copper", "Gold", "Silver", "Lead", "Tin", "Iron", "Sulfur"}

    # The "equation" is the set difference: Known Elements - Mentioned Elements = Unmentioned Elements
    # We will print each component of this operation.
    
    print("Equation: {Known Elements} - {Mentioned Elements} = {Unmentioned Elements}\n")

    print("Elements Known in Antiquity (excluding Zinc, Antimony):")
    # sorted() is used for consistent, readable output
    print(sorted(list(elements_known_in_antiquity)))
    
    print("\nElements Mentioned in the Odyssey:")
    print(sorted(list(elements_mentioned_in_odyssey)))

    # Calculate the difference between the two sets to find the unmentioned elements.
    unmentioned_elements = elements_known_in_antiquity.difference(elements_mentioned_in_odyssey)

    print("\nElements from the first list that are NOT in the second list:")
    # The loop will print each element that was not mentioned.
    for element in unmentioned_elements:
        print(element)

# Run the function to display the results.
find_unmentioned_elements()