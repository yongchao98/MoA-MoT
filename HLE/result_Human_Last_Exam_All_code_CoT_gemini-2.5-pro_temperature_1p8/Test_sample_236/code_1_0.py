import sys

def solve_homology_cobordism_question():
    """
    Calculates and explains the number of homology cobordism group elements
    representable by integral surgery on knots with at most four crossings.
    """
    # Step 1: Define the knots and their relevant properties.
    # The 'generates' key holds the set of elements in the homology cobordism group (Theta_3^H)
    # produced by +1 and -1 surgery on the knot.
    # We use descriptive strings to represent these elements.
    knot_data = {
        'Unknot (0_1)': {
            'crossings': 0,
            'is_slice': True,
            'generates': {'Identity Element (S^3)'}
        },
        'Trefoil (3_1)': {
            'crossings': 3,
            'is_slice': False,
            # +1/-1 surgeries on the (non-slice) trefoil yield the Poincaré homology sphere and its inverse.
            'generates': {'Poincaré Homology Sphere', 'Inverse Poincaré Homology Sphere'}
        },
        'Figure-eight (4_1)': {
            'crossings': 4,
            'is_slice': True,
            # Being a slice knot means +/-1 surgery results are homology cobordant to S^3.
            'generates': {'Identity Element (S^3)'}
        }
    }

    print("Analyzing integral surgeries on knots with at most 4 crossings...")
    print("A surgery on a knot K creates a homology sphere only with coefficients +1 or -1.")
    print("If a knot is 'slice', these surgeries result in the identity element in the homology cobordism group.\n")

    # Step 2: Tally the unique elements from all relevant knots.
    unique_elements = set()
    for name, properties in knot_data.items():
        print(f"Knot: {name}")
        if properties['is_slice']:
            print(f"  - Property: SLICE. Generated elements are trivial in the group.")
        else:
            print(f"  - Property: NOT SLICE. Generates non-trivial elements.")
        
        # Add the generated elements to the total set of unique elements
        unique_elements.update(properties['generates'])
        print(f"  - Contribution: {properties['generates']}\n")

    # Step 3: Summarize the findings and present the final count.
    print("---------------------------------------------------------")
    print("The complete set of unique elements generated is:")
    for i, element in enumerate(sorted(list(unique_elements))):
        print(f"  {i+1}. {element}")
    
    total_count = len(unique_elements)
    
    # Per the prompt, show the numbers that make up the final calculation.
    # We have 1 identity element and 2 non-trivial elements.
    identity_count = 1
    nontrivial_count = total_count - identity_count

    print("\nFinal Calculation:")
    print(f"   {identity_count} (the identity element, from the Unknot and Figure-eight knot)")
    print(f"+  {nontrivial_count} (non-trivial elements, from the Trefoil knot)")
    print( "--------------------")
    print(f"=  {total_count} (total distinct elements)")
    print("---------------------------------------------------------")

if __name__ == "__main__":
    solve_homology_cobordism_question()