import spherogram
import sys

# We will use the spherogram library to perform the braid analysis.
# This script will define the braid, find the components of its closure,
# and identify the knot type of each component.

def solve_braid_closure():
    """
    Analyzes the braid beta = sigma_1^2 * sigma_2^2 * sigma_3 * sigma_4^-1 in B_5
    and identifies the components of its closure.
    """
    # The braid is beta = sigma_1^2 * sigma_2^2 * sigma_3 * sigma_4^-1.
    # In spherogram notation, sigma_i is represented by the integer i,
    # and its inverse sigma_i^-1 is represented by -i.
    # The braid word is a list of these integers.
    braid_word = [1, 1, 2, 2, 3, -4]
    num_strands = 5
    
    print(f"Analyzing the closure of the braid in B_{num_strands} with word: {braid_word}")

    # Create the Braid object
    try:
        b = spherogram.Braid(word=braid_word, num_strands=num_strands)
    except ImportError as e:
        print(f"Error: {e}", file=sys.stderr)
        print("Please ensure you have installed the 'spherogram' library.", file=sys.stderr)
        print("You can install it using: pip install spherogram", file=sys.stderr)
        return

    # Compute the closure of the braid to get the corresponding link
    link = b.closure()

    # Get the list of knot components from the link
    components = link.components
    
    print(f"\nThe closure of the braid has {len(components)} components.")
    
    # Identify each component and store its type
    component_types = []
    for i, knot_component in enumerate(components):
        # The identify() method returns a list of matching knots from the database.
        # We'll take the first and most common one.
        identified_knots = knot_component.identify()
        if identified_knots:
            # The name() method gives standard knot notation, e.g., '0_1' for unknot, '3_1' for trefoil.
            knot_name = identified_knots[0].name()
            component_types.append(knot_name)
        else:
            component_types.append("Unknown")

    # Translate knot names into common English names and print the final result.
    print("\nThe components of the link are identified as follows:")
    
    unknot_count = 0
    other_component_type = "Not found"
    
    name_map = {
        '0_1': 'Unknot',
        '3_1': 'Trefoil',
        '4_1': 'Figure-8',
        '5_1': '$5_1$ knot'
    }

    for knot_id in component_types:
        knot_type = name_map.get(knot_id, f"Unknown ({knot_id})")
        print(f"- One component is a(n) {knot_type}")
        if knot_type == 'Unknot':
            unknot_count += 1
        else:
            other_component_type = knot_type
            
    print("\n--- Conclusion ---")
    if unknot_count == 2 and other_component_type != "Not found":
        print("The problem states that two components are unknots, which our analysis confirms.")
        print(f"The other connected component is equivalent to the {other_component_type} knot.")
    else:
        print("The analysis did not confirm the premise of the problem.")

if __name__ == '__main__':
    solve_braid_closure()
