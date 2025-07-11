import collections

def solve_homology_cobordism_question():
    """
    Calculates the number of homology cobordism group elements representable
    by integral surgery on knots with at most four crossings.

    The solution relies on known results from low-dimensional topology regarding
    which surgeries on specific knots produce homology spheres that are
    trivial in the homology cobordism group.
    """

    # The homology cobordism group elements are identified by names.
    # 'trivial' represents the class of the 3-sphere S^3.
    # 'poincare' represents the class of the Poincare homology sphere.
    # 'poincare_mirror' represents the class of its mirror image.
    # These three are known to be distinct elements.

    # Dictionary to store the results of +/- 1 surgery for each knot.
    # Key: knot name, Value: (result of +1 surgery, result of -1 surgery)
    surgery_results = collections.OrderedDict({
        "Unknot (0_1)": ("trivial", "trivial"),
        "Right-handed Trefoil (3_1)": ("poincare", "trivial"),
        "Left-handed Trefoil (mirror 3_1)": ("trivial", "poincare_mirror"),
        "Figure-eight Knot (4_1)": ("trivial", "trivial")
    })

    print("Step-by-step analysis of elements from surgery on knots with at most 4 crossings:")
    print("-" * 70)

    unique_elements = set()
    knot_contributions = {}

    for knot, results in surgery_results.items():
        initial_count = len(unique_elements)
        plus_one_surgery_result, minus_one_surgery_result = results
        
        # Add the elements from the surgeries to our set of unique elements
        unique_elements.add(plus_one_surgery_result)
        unique_elements.add(minus_one_surgery_result)
        
        # Calculate how many *new* elements this knot contributed
        new_elements_count = len(unique_elements) - initial_count
        knot_contributions[knot] = new_elements_count
        
        print(f"Knot: {knot}")
        print(f"  +1 surgery yields the '{plus_one_surgery_result}' element.")
        print(f"  -1 surgery yields the '{minus_one_surgery_result}' element.")
        print(f"  New elements contributed by this knot: {new_elements_count}\n")


    print("Summary of contributions:")
    total = 0
    equation_parts = []
    for knot, count in knot_contributions.items():
        print(f"  Contribution from {knot}: {count}")
        total += count
        equation_parts.append(str(count))
    
    equation = " + ".join(equation_parts)

    print("\nFinal Calculation:")
    print(f"Total number of unique elements = {equation} = {total}")

    print("\nThe distinct elements found are:")
    for element in sorted(list(unique_elements)):
        print(f"  - {element}")
    
    print(f"\nThus, the total number of distinct elements is {total}.")
    
    # Returning the final answer for the wrapper
    return total

if __name__ == "__main__":
    final_answer = solve_homology_cobordism_question()
    # The format below is for the final answer extraction.
    # print(f"<<<{final_answer}>>>")
