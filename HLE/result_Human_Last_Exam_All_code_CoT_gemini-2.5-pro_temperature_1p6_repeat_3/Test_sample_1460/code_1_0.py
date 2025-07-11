import spherogram

def solve_braid_closure():
    """
    This function analyzes the closure of the braid
    beta = sigma_1^2 * sigma_2^2 * sigma_3 * sigma_4^{-1}
    in the braid group B_5 and identifies its components.
    """
    # The braid is on n=5 strands.
    n_strands = 5
    # The braid word is represented by a list of integers.
    # sigma_i is represented by i, and sigma_i^-1 by -i.
    # beta = sigma_1^2 * sigma_2^2 * sigma_3 * sigma_4^-1 corresponds to:
    braid_word = [1, 1, 2, 2, 3, -4]

    print(f"Analyzing the braid on {n_strands} strands with word: {braid_word}")

    # Create the BraidGroup object for B_5
    B = spherogram.BraidGroup(n_strands)
    # Create the braid element from the word
    beta = B(braid_word)

    # Compute the closure of the braid to get the resulting link
    L = beta.closure()

    # The link is composed of several knot components. Let's find them.
    components = L.components
    num_components = len(components)
    print(f"The closure of the braid forms a link with {num_components} components.")

    component_names = []
    # Identify each component using spherogram's database
    for i, K in enumerate(components):
        # The identify() method returns a list of possible matches.
        # We'll take the first, most common one, e.g., Knot(4,1).
        # Its name attribute would be '4_1'. The unknot is Knot(0,1), named '0_1'.
        identity = K.identify()
        if identity:
            knot_name = identity[0].name()
            component_names.append(knot_name)
        else:
            component_names.append("Unknown")

    print(f"The identified components are: {', '.join(component_names)}")

    # According to the problem, two components are unknots ('0_1').
    unknot_count = component_names.count('0_1')
    other_components = [name for name in component_names if name != '0_1']

    print(f"\nVerifying problem statement: Found {unknot_count} unknot(s).")

    if unknot_count == 2 and len(other_components) == 1:
        final_knot_name = other_components[0]
        # Determine the common name for the knot.
        if final_knot_name == '4_1':
            final_answer = "Figure-8"
        elif final_knot_name == '3_1':
            final_answer = "Trefoil"
        elif final_knot_name == '5_1':
            final_answer = "5_1"
        elif final_knot_name == '0_1':
            final_answer = "Unknot"
        else:
            final_answer = final_knot_name

        print(f"The other connected component is the {final_knot_name} knot, also known as the {final_answer} knot.")
    else:
        print("Could not find the final answer because the conditions of the problem were not met by the analysis.")

if __name__ == "__main__":
    solve_braid_closure()