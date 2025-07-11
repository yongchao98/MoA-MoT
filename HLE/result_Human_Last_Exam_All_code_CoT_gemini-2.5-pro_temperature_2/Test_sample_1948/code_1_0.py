def find_unmentioned_elements():
    """
    This script identifies which anciently-known chemical elements
    are not mentioned in Homer's The Odyssey.
    """

    # The set of elements known in antiquity, excluding zinc and antimony as per the prompt.
    ancient_elements = {
        "Gold",
        "Silver",
        "Copper",
        "Iron",
        "Lead",
        "Tin",
        "Mercury",
        "Carbon",
        "Sulfur"
    }

    # The set of the above elements that are mentioned in The Odyssey.
    # Mentions can be direct or through common forms like bronze (Copper), soot (Carbon), etc.
    mentioned_in_odyssey = {
        "Gold",
        "Silver",
        "Copper",
        "Iron",
        "Lead",
        "Tin",
        "Carbon",
        "Sulfur"
    }

    # The "equation" is finding the difference between the two sets.
    # We will output all the elements involved in this logical operation.
    unmentioned_elements = ancient_elements.difference(mentioned_in_odyssey)

    print("Step 1: Define the set of all relevant ancient elements.")
    print("All Elements = ", sorted(list(ancient_elements)))
    print("-" * 20)

    print("Step 2: Define the set of elements mentioned in The Odyssey.")
    print("Mentioned Elements = ", sorted(list(mentioned_in_odyssey)))
    print("-" * 20)

    print("Step 3: Calculate the difference to find the unmentioned elements.")
    print("Result (All Elements - Mentioned Elements) = ", list(unmentioned_elements))
    print("-" * 20)
    
    if unmentioned_elements:
        print("\nThe element from the list of anciently known elements that is NOT mentioned in The Odyssey is:")
        for element in unmentioned_elements:
            print(element)
    else:
        print("All ancient elements were mentioned in The Odyssey.")

find_unmentioned_elements()