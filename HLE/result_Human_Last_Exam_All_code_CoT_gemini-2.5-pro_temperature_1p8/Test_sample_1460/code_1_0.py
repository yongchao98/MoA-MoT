import spherogram as sg

def solve_braid_closure():
    """
    This function defines a braid, computes its closure,
    and identifies the knot type of each resulting component.
    """
    try:
        # 1. Represent the braid from the problem statement.
        # The braid is beta = sigma_1^2 * sigma_2^2 * sigma_3 * sigma_4^-1 in B_5.
        # In spherogram's notation, sigma_i is represented by the integer i,
        # and its inverse by -i.
        braid_word = [1, 1, 2, 2, 3, -4]
        braid = sg.Braid(5, braid_word)
        print(f"Analyzing the closure of the braid: {braid}")

        # 2. Compute the closure of the braid to get a link.
        link = braid.closure()

        # 3. Decompose the link into its components.
        components = link.components
        print(f"\nThe resulting link has {len(components)} components.")
        print("Identifying the knot type of each component...")

        # 4. Identify each component's knot type.
        component_types = []
        for i, c in enumerate(components):
            # The identify() method returns a list of possible knots,
            # with the most likely first.
            identified_knot = c.identify()[0]
            component_types.append(identified_knot)
            print(f"  Component {i+1} is a(n) {identified_knot.name()}.")
        
        # 5. Find the non-unknot component and match it with the options.
        unknot_count = 0
        final_knot = None
        for knot in component_types:
            # The unknot is named '0_1'.
            if knot.name() == '0_1':
                unknot_count += 1
            else:
                final_knot = knot

        print(f"\nAs stated in the problem, we found {unknot_count} unknots.")

        if final_knot:
            knot_name = final_knot.name()
            answer_name = ""
            if knot_name == '3_1':
                answer_name = "Trefoil"
            elif knot_name == '4_1':
                answer_name = "Figure-8"
            elif knot_name == '5_1':
                answer_name = "5_1 knot"
            
            if answer_name:
                 print(f"The other component is the {knot_name} knot, also known as the {answer_name} knot.")
            else:
                 print(f"The other component is the {knot_name} knot.")
        else:
            print("No non-unknot component was found.")

    except ImportError:
        print("This script requires the 'spherogram' library.")
        print("Please install it using: pip install spherogram")
    except Exception as e:
        print(f"An error occurred: {e}")

solve_braid_closure()