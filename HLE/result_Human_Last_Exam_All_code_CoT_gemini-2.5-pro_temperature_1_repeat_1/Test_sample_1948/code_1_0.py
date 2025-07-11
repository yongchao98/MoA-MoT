def find_unmentioned_elements():
    """
    Identifies which anciently known chemical elements are not mentioned in the Odyssey.
    """
    # Elements known in elemental form around the time of the Odyssey's composition
    # (excluding zinc and antimony as per the user's request).
    ancient_elements = {"Gold", "Silver", "Copper", "Tin", "Lead", "Iron", "Mercury", "Carbon", "Sulfur"}

    # Elements mentioned in the Odyssey, either directly or as a primary component of a material.
    # 'Bronze' implies copper and tin. 'Soot' implies carbon. 'Brimstone' is sulfur.
    mentioned_in_odyssey = {"Gold", "Silver", "Copper", "Tin", "Lead", "Iron", "Carbon", "Sulfur"}

    # Find the difference between the two sets to identify the unmentioned elements.
    unmentioned_elements = ancient_elements.difference(mentioned_in_odyssey)

    print("Elements known in antiquity (excluding Zn, Sb):")
    print(sorted(list(ancient_elements)))
    print("\nElements mentioned in the Odyssey:")
    print(sorted(list(mentioned_in_odyssey)))
    print("\nElements from the first list that are NOT mentioned in the Odyssey:")
    for element in unmentioned_elements:
        print(element)

find_unmentioned_elements()
<<<Mercury>>>