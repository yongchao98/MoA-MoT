def find_bacterium_by_property(phenomenon):
    """
    Simulates a query to a biological knowledge base to find a bacterium
    genus based on a specific chemical property.
    """
    # A small knowledge base mapping phenomena to bacterial genera.
    knowledge_base = {
        'triboluminescence': 'Photorhabdus',
        'bioluminescence': 'Vibrio',
        'nitrogen_fixation': 'Rhizobium',
        'spore_formation': 'Bacillus'
    }

    # Retrieve the genus associated with the given phenomenon.
    genus = knowledge_base.get(phenomenon, "not found in this knowledge base")

    # Print the result in a clear format.
    print(f"Searching for a bacterium with the property: '{phenomenon}'")
    if genus != "not found in this knowledge base":
        print(f"The genus of bacterium that produces a chemical compound exhibiting {phenomenon} properties is: {genus}")
    else:
        print("Could not find a bacterium with that property.")

# The specific property we are interested in.
search_property = 'triboluminescence'
find_bacterium_by_property(search_property)