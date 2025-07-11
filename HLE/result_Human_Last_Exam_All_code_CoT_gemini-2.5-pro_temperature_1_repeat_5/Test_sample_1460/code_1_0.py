import sys
# It's good practice to check for pyknotid and provide a helpful message if it's not found.
try:
    from pyknotid.spacecurves import Braid
    from pyknotid.catalogue import identify
except ImportError:
    print("Please install the 'pyknotid' library to run this script.")
    print("You can install it using: pip install pyknotid")
    sys.exit(1)

def solve_braid_closure():
    """
    This function solves the problem by identifying the components of the closure
    of the braid beta = sigma_1^2 * sigma_2^2 * sigma_3 * sigma_4^-1 in B_5.
    """
    # 1. Define the braid word for beta in B_5.
    # B_5 has 5 strands and 4 generators (sigma_1 to sigma_4).
    # The braid word sigma_1^2 * sigma_2^2 * sigma_3 * sigma_4^-1 is represented as a list of integers.
    braid_word = [1, 1, 2, 2, 3, -4]
    num_strands = 5

    print(f"Analyzing the closure of the braid on {num_strands} strands with word: {braid_word}")

    # 2. Create the Braid object and compute its closure.
    try:
        braid = Braid(braid_word, num_strands=num_strands)
        link_closure = braid.closure()
    except Exception as e:
        print(f"An error occurred during braid closure computation: {e}")
        return

    # 3. Get the individual components of the resulting link.
    components = link_closure.components()
    num_components = len(components)
    print(f"The closure resulted in a link with {num_components} components.")

    if num_components != 3:
        print("Warning: The number of components is not 3 as implied by the problem statement.")

    # 4. Identify the knot type of each component.
    component_names = []
    component_map = {
        '0_1': 'Unknot',
        '3_1': 'Trefoil',
        '4_1': 'Figure-8',
        '5_1': '5_1 knot'
    }

    print("\n--- Identifying Components ---")
    for i, component in enumerate(components):
        # The identify() function returns a list of possible matches.
        # We take the first and most likely one.
        identity = identify(component)
        if identity:
            # The result is a tuple (Knot, periodic, mirrored)
            knot_id_str = str(identity[0][0])
            component_names.append(knot_id_str)
            common_name = component_map.get(knot_id_str, knot_id_str) # Use common name if available
            print(f"Component {i + 1} is a(n) {common_name} (Rolfsen notation: {knot_id_str}).")
        else:
            component_names.append('Unknown')
            print(f"Component {i + 1} could not be identified.")
    
    # 5. Determine the answer based on the problem statement.
    print("\n--- Conclusion ---")
    unknot_id = '0_1'
    num_unknots = component_names.count(unknot_id)
    
    if num_unknots == 2:
        print("As stated in the problem, two components are unknots.")
        other_components = [name for name in component_names if name != unknot_id]
        if len(other_components) == 1:
            other_knot_id = other_components[0]
            other_knot_name = component_map.get(other_knot_id, other_knot_id)
            print(f"The other connected component is the {other_knot_name}.")
        else:
            print("Error: Expected exactly one non-unknot component, but found a different number.")
    else:
        print(f"Warning: Found {num_unknots} unknots, which contradicts the problem statement.")


if __name__ == "__main__":
    solve_braid_closure()