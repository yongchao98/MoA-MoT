def find_unstable_complexes():
    """
    Identifies iridium complexes expected to have shorter lifetimes based on their structure.

    The key factor for instability in these complexes is the presence of a fluorine atom
    at the 2'-position of the phenylpyridine ligand (ortho to the Ir-C bond),
    which introduces a known degradation pathway.
    """

    # Data representing the four complexes and whether they have the destabilizing 2'-fluoro group.
    # True means the feature is present, leading to a shorter expected lifetime.
    complex_data = {
        1: {'has_destabilizing_feature': False, 'description': "Fluorines at 4' and 5' positions"},
        2: {'has_destabilizing_feature': False, 'description': "No fluorine substituents"},
        3: {'has_destabilizing_feature': True,  'description': "Fluorine at 2' position"},
        4: {'has_destabilizing_feature': True,  'description': "Fluorines at 2' and 4' positions"}
    }

    unstable_complex_numbers = []
    for complex_id, data in complex_data.items():
        if data['has_destabilizing_feature']:
            unstable_complex_numbers.append(complex_id)

    print("Based on chemical principles, complexes with a fluorine atom at the 2'-position of the phenyl ring are prone to degradation and exhibit shorter lifetimes.")
    print("The complexes that possess this structural feature are:")
    
    # Print each number as requested
    output = " and ".join([f"Complex {num}" for num in unstable_complex_numbers])
    print(output)
    
    # To be explicit with the numbers for the final answer
    print("\nFinal identified numbers:")
    for num in unstable_complex_numbers:
        print(num)

find_unstable_complexes()