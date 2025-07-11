def find_unmentioned_elements():
    """
    Identifies which of the anciently-known chemical elements are not mentioned in the Odyssey.
    """
    # Elements known in their elemental form in antiquity (c. 8th century BC),
    # excluding zinc and antimony as per the user's request.
    elements_known_in_antiquity = {"Carbon", "Sulfur", "Iron", "Copper", "Silver", "Tin", "Gold", "Mercury", "Lead"}

    # Elements from the above list that are mentioned in the Odyssey.
    # Copper and Tin are mentioned as components of bronze.
    # Carbon is mentioned as soot from the fire.
    # Sulfur is mentioned by name for purifying the great hall.
    elements_mentioned_in_odyssey = {"Gold", "Silver", "Copper", "Tin", "Iron", "Lead", "Sulfur", "Carbon"}

    # Find the difference between the two sets to identify the unmentioned elements.
    unmentioned_elements = elements_known_in_antiquity.difference(elements_mentioned_in_odyssey)

    # Print the result.
    if unmentioned_elements:
        print("The element known to the ancient world but not mentioned in the Odyssey is:")
        for element in unmentioned_elements:
            print(element)
    else:
        print("All elements known in antiquity were mentioned in the Odyssey.")

find_unmentioned_elements()