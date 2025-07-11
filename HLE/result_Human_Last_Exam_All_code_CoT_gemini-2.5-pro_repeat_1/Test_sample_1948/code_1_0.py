def find_unmentioned_elements():
    """
    Finds which anciently known chemical elements are not mentioned in Homer's Odyssey.
    """
    # Step 1: Define the set of elements known in elemental form around the time of the Odyssey's composition.
    # Zinc and Antimony are excluded as per the prompt.
    ancient_elements = {
        'Carbon', 'Copper', 'Gold', 'Silver', 'Lead', 'Tin', 'Iron', 'Sulfur', 'Mercury'
    }

    # Step 2: Define the set of elements from the above list that are mentioned in the Odyssey.
    # Gold (Au), Silver (Ag), Copper (Cu), Tin (Sn), Iron (Fe), Lead (Pb), Sulfur (S), Carbon (C)
    odyssey_elements = {
        'Gold', 'Silver', 'Copper', 'Tin', 'Iron', 'Lead', 'Sulfur', 'Carbon'
    }

    # Step 3: Find the difference between the two sets to identify the unmentioned elements.
    unmentioned_elements = ancient_elements.difference(odyssey_elements)

    # Step 4: Print the result.
    if unmentioned_elements:
        print("Of the elements known in their elemental form at the time (excluding zinc and antimony), the following are not mentioned in the Odyssey:")
        for element in unmentioned_elements:
            print(f"- {element}")
    else:
        print("All anciently known elements (excluding zinc and antimony) are mentioned in the Odyssey.")

find_unmentioned_elements()
<<<Mercury>>>