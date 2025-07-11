import spherogram as sp

def solve_braid_closure():
    """
    This function solves the problem by constructing the braid in spherogram,
    finding its closure, and identifying the knot type of each component.
    """
    try:
        # 1. Define the braid group B_5 and the braid beta.
        # In spherogram, generators are indexed from 0, so sigma_i corresponds to gen(i-1).
        B5 = sp.BraidGroup(5)
        s = B5.generators()
        # beta = sigma_1^2 * sigma_2^2 * sigma_3 * sigma_4^{-1}
        beta = s[0]**2 * s[1]**2 * s[2] * s[3]**-1

        # This is not a real equation, but to satisfy the prompt's instruction to
        # "output each number in the final equation", we print the powers in the braid word.
        print("The braid is beta = sigma_1^(2) * sigma_2^(2) * sigma_3^(1) * sigma_4^(-1)")
        print("Powers in the braid expression: 2, 2, 1, -1")
        print("-" * 30)

        # 2. Compute the link from the closure of the braid.
        L = sp.Link(beta)

        # 3. Get the components of the link.
        components = L.components
        num_components = len(components)
        print(f"The closure of the braid forms a link with {num_components} components.")

        # 4. Identify the knot type of each component.
        component_knots = []
        for i, C in enumerate(components):
            # To identify a component, we create a new Link object from it.
            knot = sp.Link(C)
            knot_name = knot.name()

            # Convert standard names to common names for clarity.
            if knot_name == '0_1':
                common_name = 'Unknot'
            elif knot_name == '3_1':
                common_name = 'Trefoil'
            elif knot_name == '4_1':
                common_name = 'Figure-8'
            elif knot_name == '5_1':
                common_name = '5_1 knot'
            else:
                common_name = f"Unknown ({knot_name})"
            component_knots.append(common_name)
            print(f"Component {i+1} is a(n): {common_name}")

        # 5. Determine the answer based on the problem statement.
        other_component = "Not Found"
        for knot_type in component_knots:
            if knot_type != 'Unknot':
                other_component = knot_type
                break
        
        print("-" * 30)
        print("Given that two of the components are unknots, the other connected component is the:")
        print(other_component)

    except ImportError:
        print("This script requires the 'spherogram' library.")
        print("Please install it using: pip install spherogram")
    except Exception as e:
        print(f"An error occurred: {e}")
        print("This may be due to a missing dependency like 'snappy'.")
        print("Please ensure spherogram and its dependencies are correctly installed.")

if __name__ == "__main__":
    solve_braid_closure()