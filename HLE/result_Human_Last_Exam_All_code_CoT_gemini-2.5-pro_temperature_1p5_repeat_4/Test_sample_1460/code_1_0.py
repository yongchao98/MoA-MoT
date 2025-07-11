import sys
try:
    import spherogram
except ImportError:
    # If spherogram is not installed, print an installation message and exit.
    print("This script requires the 'spherogram' library.")
    print("Please install it, for example using pip: 'pip install spherogram'")
    sys.exit(1)

def solve_knot_problem():
    """
    This function solves the given braid problem by:
    1. Defining the braid in the braid group B_5.
    2. Taking its closure to form a link.
    3. Identifying the components of the link.
    4. Determining the knot type of the non-unknot component.
    """

    # Define the braid group on 5 strands, B_5
    B5 = spherogram.BraidGroup(5)
    
    # Get the standard generators for B_5
    s1, s2, s3, s4 = B5.generators()

    # The braid is beta = sigma_1^2 * sigma_2^2 * sigma_3 * sigma_4^-1
    beta = s1**2 * s2**2 * s3 * s4**-1

    print(f"Analyzing the closure of the braid: {beta}")

    # Take the closure of the braid to get the link L
    L = beta.closure()

    # Get the connected components of the link
    components = L.link_components
    
    print(f"The resulting link has {len(components)} components.")
    print("Given that two components are unknots, we will identify the third.")
    
    # Find the component that is not an unknot
    non_unknot_component = None
    for K in components:
        if not K.is_unknot():
            non_unknot_component = K
            break
            
    # Check if a non-unknot component was found
    if non_unknot_component:
        # Use the .name() method to get the standard knot identifier (Rolfsen notation)
        knot_name = non_unknot_component.name()
        
        # Map the standard name to the common name
        if knot_name == '4_1':
            common_name = "Figure-8"
        elif knot_name == '3_1':
            common_name = "Trefoil"
        elif knot_name == '5_1':
            common_name = "5_1 knot"
        else:
            common_name = f"knot with identifier {knot_name}"
            
        print(f"\nThe third component is the {common_name} knot.")
    else:
        print("\nAll components were identified as unknots.")

if __name__ == "__main__":
    solve_knot_problem()
