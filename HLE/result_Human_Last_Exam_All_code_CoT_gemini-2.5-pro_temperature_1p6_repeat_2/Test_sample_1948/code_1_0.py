def find_unmentioned_elements():
    """
    This script identifies which chemical elements known in antiquity (at the time of
    The Odyssey's composition) are not mentioned in the poem.
    """
    
    # Step 1: Define the set of elements known in antiquity, excluding zinc and antimony.
    ancient_elements = {'Carbon', 'Copper', 'Gold', 'Iron', 'Lead', 'Mercury', 'Silver', 'Sulfur', 'Tin'}
    
    # Step 2: Define the set of the above elements that are mentioned in The Odyssey.
    # Copper is mentioned as bronze, and Carbon as charcoal.
    odyssey_elements = {'Carbon', 'Copper', 'Gold', 'Iron', 'Lead', 'Silver', 'Sulfur', 'Tin'}
    
    # Step 3: Calculate the difference to find the elements not mentioned.
    not_mentioned_elements = ancient_elements.difference(odyssey_elements)
    
    # For clear output, convert sets to sorted lists.
    ancient_list_str = " + ".join(sorted(list(ancient_elements)))
    odyssey_list_str = " + ".join(sorted(list(odyssey_elements)))
    result_list_str = " + ".join(sorted(list(not_mentioned_elements)))

    print("This script solves which anciently-known elements are missing from The Odyssey.")
    print("-" * 70)
    print("Elements known in antiquity:")
    print(f"[{ancient_list_str}]")
    print("\nElements mentioned in The Odyssey:")
    print(f"[{odyssey_list_str}]")
    print("-" * 70)

    # Step 4: Print the final result as a clear equation, as requested.
    print("Final Equation:")
    print(f"({ancient_list_str}) - ({odyssey_list_str}) = {result_list_str}")
    print("-" * 70)

    print("\nThe element known in antiquity but NOT mentioned in The Odyssey is:")
    for element in not_mentioned_elements:
        print(element)

find_unmentioned_elements()
<<<Mercury>>>