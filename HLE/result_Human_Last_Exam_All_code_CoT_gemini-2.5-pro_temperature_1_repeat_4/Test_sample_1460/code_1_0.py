import spherogram as sp

def solve_knot_problem():
    """
    This script identifies the components of a link formed by the closure of a specific braid.
    The braid is beta = sigma_1^2 * sigma_2^2 * sigma_3 * sigma_4^-1 in B_5.
    We identify the components by computing and comparing their Jones polynomials.
    """
    # 1. Define the braid from the problem statement.
    # In spherogram, a braid on n strands is defined by a list of integers.
    # A positive integer i corresponds to the generator sigma_i.
    # A negative integer -i corresponds to the inverse generator sigma_i^-1.
    num_strands = 5
    # The braid beta = sigma_1^2 * sigma_2^2 * sigma_3 * sigma_4^-1 is represented as:
    braid_list = [1, 1, 2, 2, 3, -4]
    
    print(f"Analyzing the link from the closure of the braid on {num_strands} strands.")
    print(f"The braid's generator sequence is: {braid_list}\n")

    # 2. Create the braid object and compute its closure to get the link.
    braid = sp.Braid(num_strands, braid_list)
    link = sp.Link(braid)

    # 3. Decompose the link into its connected components.
    components = link.split()
    print(f"The resulting link has {len(components)} components.")
    print("According to the problem, two components are unknots. Let's find the third.\n")

    # 4. Calculate and print the Jones polynomial for each component.
    print("--- Component Analysis ---")
    component_polynomials = []
    for i, c in enumerate(components):
        jp = c.jones_polynomial()
        component_polynomials.append(jp)
        print(f"Component {i+1} has Jones Polynomial: {jp}")

    print("\n--- Reference Knot Polynomials ---")
    # 5. Compute the Jones polynomials for the reference knots from the answer choices.
    ref_knots = {
        "Unknot": sp.Link('0_1'),
        "Trefoil": sp.Link('3_1'),
        "Figure-8": sp.Link('4_1'),
        "$5_1$ knot": sp.Link('5_1')
    }

    ref_polynomials = {}
    for name, knot in ref_knots.items():
        jp = knot.jones_polynomial()
        ref_polynomials[name] = jp
        print(f"{name}: {jp}")

    # 6. Identify the non-unknot component by comparing polynomials.
    print("\n--- Conclusion ---")
    unknot_poly = ref_polynomials["Unknot"]
    non_unknot_components = []
    for poly in component_polynomials:
        if poly != unknot_poly:
            non_unknot_components.append(poly)

    if len(non_unknot_components) == 1:
        the_knot_poly = non_unknot_components[0]
        found = False
        for name, poly in ref_polynomials.items():
            if the_knot_poly == poly:
                print(f"The non-unknot component has a polynomial that matches the {name}.")
                found = True
                break
        if not found:
            print("The non-unknot component does not match the provided choices.")
    else:
        print("The link does not have two unknot components as stated.")

if __name__ == '__main__':
    # You may need to install the spherogram library first:
    # pip install spherogram
    solve_knot_problem()