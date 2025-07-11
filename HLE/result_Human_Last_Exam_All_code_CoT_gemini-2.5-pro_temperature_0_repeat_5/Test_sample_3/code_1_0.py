import sys

def solve_torsion_problem():
    """
    Finds the largest order of a non-cyclic torsion subgroup of an elliptic curve
    over the number field Q(sqrt(-3)) by applying the classification theorem for
    torsion subgroups over quadratic fields.
    """
    # Step 1 & 2: The complete list of possible non-cyclic torsion subgroups over
    # quadratic fields, based on the classification theorem.
    # The data is stored as a list of tuples: (Group Name, Order, Required Discriminant(s)).
    # 'any' means the group can appear over Q, and thus any quadratic field.
    # Otherwise, the list contains the specific field discriminant(s) required.
    # The field Q(sqrt(-3)) has discriminant -3.
    known_groups = [
        ("Z/2Z x Z/2Z", 4, 'any'),
        ("Z/2Z x Z/4Z", 8, 'any'),
        ("Z/2Z x Z/6Z", 12, 'any'),
        ("Z/2Z x Z/8Z", 16, 'any'),
        ("Z/3Z x Z/3Z", 9, [-3]),
        ("Z/3Z x Z/6Z", 18, [-3]),
        ("Z/4Z x Z/4Z", 16, [-4]),
        ("Z/2Z x Z/10Z", 20, [5]),
        ("Z/2Z x Z/12Z", 24, [-4]),
    ]

    target_discriminant = -3
    
    print(f"The problem is to find the largest order of a non-cyclic torsion subgroup of an elliptic curve over Q(sqrt(-3)).")
    print(f"This field has a discriminant of {target_discriminant}.")
    print("We use the known classification of torsion subgroups over quadratic fields.")
    print("\nThe possible non-cyclic torsion groups over Q(sqrt(-3)) are:")

    # Step 3: Filter the list to find groups possible over Q(sqrt(-3)).
    possible_orders = []
    for name, order, d_list in known_groups:
        if d_list == 'any' or target_discriminant in d_list:
            possible_orders.append(order)
            print(f"- The group {name}, which has order {order}.")

    # Step 4: Find the maximum order from the filtered list.
    if not possible_orders:
        # This case should not happen based on the known theory.
        max_order = 0
        equation_str = "No possible non-cyclic groups found."
    else:
        max_order = max(possible_orders)
        # Create a string representing the final calculation as requested.
        sorted_orders = sorted(list(set(possible_orders)))
        equation_str = f"max({', '.join(map(str, sorted_orders))}) = {max_order}"

    print("\nTo find the largest order, we take the maximum of the possible orders:")
    print(equation_str)

if __name__ == "__main__":
    solve_torsion_problem()
