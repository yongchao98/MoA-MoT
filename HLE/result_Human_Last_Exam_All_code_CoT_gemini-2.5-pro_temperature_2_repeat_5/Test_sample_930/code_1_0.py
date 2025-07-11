def find_milkweed_mutualists():
    """
    Identifies the mutualists of Asclepias fascicularis from a predefined list.

    The function uses a hardcoded list of organisms and their relationships
    with the milkweed plant to determine which ones are mutualists.
    Mutualism is defined as a relationship beneficial to both organisms,
    such as pollination or protection.
    """
    
    # A list of dictionaries representing the organisms and their ecological role.
    # Relationships are categorized as 'mutualism', 'herbivory', or 'none'.
    organisms = [
        {'id': 1, 'name': 'Danaus plexipus (Adult)', 'relationship': 'mutualism'},
        {'id': 2, 'name': 'Megachile frigidus (Adult)', 'relationship': 'mutualism'},
        {'id': 3, 'name': 'Formica rufa (Adult)', 'relationship': 'mutualism'},
        {'id': 4, 'name': 'Sphex ichneumoneus (Adult)', 'relationship': 'mutualism'},
        {'id': 5, 'name': 'Pepsis thisbe (Adult)', 'relationship': 'mutualism'},
        {'id': 6, 'name': 'Megachile ericetorum (Adult)', 'relationship': 'mutualism'},
        {'id': 7, 'name': 'Danaus plexipus (Larva)', 'relationship': 'herbivory'},
        {'id': 8, 'name': 'Megachile frigidus (Larva)', 'relationship': 'none'},
        {'id': 9, 'name': 'Formica rufa (Larva)', 'relationship': 'none'},
        {'id': 10, 'name': 'Sphex ichneumoneus (Larva)', 'relationship': 'none'},
        {'id': 11, 'name': 'Pepsis thisbe (Larva)', 'relationship': 'none'},
        {'id': 12, 'name': 'Megachile ericetorum (Larva)', 'relationship': 'none'}
    ]

    mutualist_indices = []
    for org in organisms:
        if org['relationship'] == 'mutualism':
            # Append the ID as a string for joining later
            mutualist_indices.append(str(org['id']))
            
    if mutualist_indices:
        result = ",".join(mutualist_indices)
    else:
        result = "none"
        
    print(result)

# Execute the function to find and print the result.
find_milkweed_mutualists()