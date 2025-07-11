import collections

def solve_homology_cobordism_elements():
    """
    Calculates the number of homology cobordism group elements from knots
    with at most four crossings.
    """
    # Knot data: A list of tuples, where each tuple contains:
    # (Knot Name, Arf Invariant, a2 Invariant)
    # a2 is the second coefficient of the Alexander-Conway polynomial, Delta''(1)/2.
    knots = [
        ("Unknot (0_1)", 0, 0),
        ("Trefoil knot (3_1)", 1, 1),
        ("Figure-eight knot (4_1)", 0, 1),
    ]

    # Use a set to store unique element descriptions to avoid duplicates.
    unique_elements = set()
    element_counts = collections.defaultdict(int)

    print("Analyzing integral surgeries on knots up to 4 crossings:")
    print("-" * 60)

    for name, arf, a2 in knots:
        print(f"Processing {name}: Arf={arf}, a2={a2}")
        if arf == 1:
            # An Arf invariant of 1 corresponds to the unique order-2 element
            # in the homology cobordism group (the Poincaré homology sphere).
            # Both +1 and -1 surgeries on the trefoil (and its mirror) yield this element.
            element_desc = "The non-trivial element of order 2 (from Trefoil knot)"
            if element_desc not in unique_elements:
                unique_elements.add(element_desc)
                element_counts["trefoil"] = 1
            print(" -> Result: A non-trivial element of order 2 (Rokhlin invariant 1).\n")
        elif arf == 0:
            if a2 == 0:
                # Arf=0 and a2=0 means the surgery results in a manifold homology
                # cobordant to S^3, the trivial element.
                element_desc = "The trivial element (from Unknot)"
                if element_desc not in unique_elements:
                    unique_elements.add(element_desc)
                    element_counts["unknot"] = 1
                print(" -> Result: The trivial element.\n")
            else: # a2 != 0
                # Arf=0 and a2!=0 means the surgery results in non-trivial elements
                # of infinite order. The +1 and -1 surgeries give an element
                # and its distinct inverse.
                element_desc_1 = "A non-trivial element of infinite order (from Figure-eight knot)"
                element_desc_2 = "The inverse of the above element (from Figure-eight knot)"
                if element_desc_1 not in unique_elements:
                    unique_elements.add(element_desc_1)
                    unique_elements.add(element_desc_2)
                    element_counts["figure-eight"] = 2
                print(" -> Result: Two distinct, non-trivial elements of infinite order.\n")

    print("-" * 60)
    print("Summary of distinct elements found:")
    for i, elem in enumerate(sorted(list(unique_elements))):
        print(f"{i+1}. {elem}")

    count_unknot = element_counts["unknot"]
    count_trefoil = element_counts["trefoil"]
    count_figure_eight = element_counts["figure-eight"]
    total_count = count_unknot + count_trefoil + count_figure_eight

    print("\nFinal Calculation:")
    print(f"Total Elements = (from Unknot) + (from Trefoil) + (from Figure-eight)")
    print(f"{total_count} = {count_unknot} + {count_trefoil} + {count_figure_eight}")


solve_homology_cobordism_elements()
<<<4>>>