def find_unmentioned_elements():
    """
    This function identifies chemical elements known in antiquity but not mentioned in the Odyssey.
    """
    # Elements known in elemental form around the 8th century BC, excluding zinc and antimony.
    known_in_antiquity = {"Gold", "Silver", "Copper", "Tin", "Lead", "Iron", "Mercury", "Sulfur", "Carbon"}

    # Elements mentioned in the Odyssey, either by name or as a common material.
    # Bronze is an alloy of Copper and Tin. Charcoal is a form of Carbon.
    mentioned_in_odyssey = {"Gold", "Silver", "Iron", "Sulfur", "Lead", "Copper", "Tin", "Carbon"}

    # Find the elements in the first set that are not in the second set.
    unmentioned_elements = known_in_antiquity.difference(mentioned_in_odyssey)

    print("Elements known to the ancients (excluding Zinc and Antimony):")
    print(", ".join(sorted(list(known_in_antiquity))))
    print("\nElements mentioned in the Odyssey:")
    print(", ".join(sorted(list(mentioned_in_odyssey))))
    
    if unmentioned_elements:
        print("\nOf the elements known in their elemental form at the time the Odyssey was composed, the element not mentioned in the poem is:")
        # The loop ensures we print the final answer directly as requested.
        for element in unmentioned_elements:
            print(f">>>{element}<<<")
    else:
        print("\nAll known elements from that period (excluding Zinc and Antimony) were mentioned in the Odyssey.")

find_unmentioned_elements()