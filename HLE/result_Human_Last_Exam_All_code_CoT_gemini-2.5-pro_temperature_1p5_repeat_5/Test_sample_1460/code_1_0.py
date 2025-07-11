import spherogram

def solve_braid_closure():
    """
    This function solves the problem by representing the given braid,
    computing its closure, and identifying the resulting link components.
    """
    # 1. Define the braid beta in B_5
    # The braid is beta = sigma_1^2 * sigma_2^2 * sigma_3 * sigma_4^{-1}.
    # In spherogram notation, sigma_i is represented by i and sigma_i^{-1} by -i.
    braid_word = [1, 1, 2, 2, 3, -4]
    try:
        braid = spherogram.Braid(braid_word)
    except ImportError as e:
        print("Spherogram library not found. Please install it using 'pip install spherogram'")
        print("Additionally, you might need to install snappy: 'pip install snappy'")
        return

    print(f"Analyzing the braid beta defined by the word: {braid_word}")
    
    # 2. Compute the closure of the braid to get a link.
    link = braid.closure()

    # 3. Get the components of the link.
    # The number of components is determined by the number of cycles
    # in the braid's permutation.
    components = link.components
    print(f"\nThe closure of the braid results in a link with {len(components)} components.")

    # 4. Identify the knot type of each component.
    print("The components are identified as:")
    knot_names = []
    for i, knot_component in enumerate(components):
        # identify() returns a list of possible matches. We take the first one.
        identified_knots = knot_component.identify()
        if identified_knots:
            # The result is a list of Knot objects, so we get their names.
            name = identified_knots[0].name()
        else:
            # An empty list from identify() usually means it's the unknot.
            # We can verify with is_unknot().
            if knot_component.is_unknot():
                name = "Unknot (0_1)"
            else:
                # This case is unlikely for simple knots.
                name = "Unknown"
        knot_names.append(name)
        print(f"  - Component {i+1}: {name}")

    # 5. Determine the final answer.
    # The problem states two components are unknots. We find the other one.
    non_unknot_components = [name for name in knot_names if "0_1" not in name]

    print("\n" + "="*40)
    print("Final Conclusion:")
    
    if len(non_unknot_components) == 1:
        final_knot = non_unknot_components[0]
        if "3_1" in final_knot:
            answer = "Trefoil"
        elif "4_1" in final_knot:
            answer = "Figure-8"
        elif "5_1" in final_knot:
            answer = "$5_1$ knot"
        else:
            answer = f"a knot known as {final_knot}"
            
        print(f"The link is composed of two Unknots and one {final_knot} knot.")
        print(f"The third connected component is equivalent to the {answer} knot.")
        
    elif len(non_unknot_components) == 0:
        answer = "Unknot"
        print("All three components are unknots.")
        print(f"The third connected component is also an {answer}.")
    else:
        print("The analysis found more than one non-trivial knot, which contradicts the problem statement.")
        answer = "Inconclusive"

    # Final "equation" output format
    print("\nThe components of the closure of B = sigma_1^2 * sigma_2^2 * sigma_3 * sigma_4^-1 are:")
    for i, name in enumerate(knot_names):
        display_name = name.replace(" (0_1)", "").replace(" (3_1)", "") # Clean up names for final display
        print(f"  Component {i+1} is the {display_name} knot.")
    
    final_answer = non_unknot_components[0].split(' ')[0] if non_unknot_components else "Unknot"
    print(f"\nThe non-unknot component is: {final_answer}")
    
solve_braid_closure()