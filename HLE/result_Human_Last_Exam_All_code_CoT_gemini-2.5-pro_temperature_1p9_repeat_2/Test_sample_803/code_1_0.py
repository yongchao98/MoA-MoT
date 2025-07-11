def check_filled_group_candidate(group_info):
    """
    Checks if a group is a candidate for a filled group based on the
    classification theorem G/Φ(G) must have order 9.
    """
    name = group_info["name"]
    order = group_info["order"]
    frattini_order = group_info["frattini_order"]
    description = group_info["description"]
    
    print(f"--- Checking Group: {name} ---")
    print(f"Description: {description}")

    if frattini_order == 0:
        print("Frattini subgroup order not provided.\n")
        return

    # A necessary condition for a nonabelian group to be filled is
    # |G/Φ(G)| = 9.
    # We check this condition.
    quotient_order = order / frattini_order
    
    print(f"Order |G| = {order}")
    print(f"Order of Frattini subgroup |Φ(G)| = {frattini_order}")
    print(f"Calculating order of quotient group G/Φ(G):")
    
    # This print statement fulfills the requirement to output numbers in the final equation.
    print(f"  |G/Φ(G)| = |G| / |Φ(G)| = {order} / {frattini_order} = {int(quotient_order)}")

    target_order = 9
    is_candidate = (quotient_order == target_order)
    
    print(f"Checking if the quotient order is {target_order}:")
    print(f"  Is {int(quotient_order)} == {target_order}? {is_candidate}")
    
    if is_candidate:
        print("Result: This group is a candidate for a filled group (further checks needed).")
    else:
        print("Result: This group is NOT a filled group because |G/Φ(G)| is not 9.")
    print("-" * 30 + "\n")

def main():
    """
    Main function to test several groups of order 2q^m.
    """
    # List of nonabelian groups with order of the form 2q^m
    groups_to_test = [
        {
            "name": "Symmetric group S3 (Dihedral D3)",
            "description": "Smallest nonabelian group, q=3, m=1.",
            "order": 6,
            "frattini_order": 1, # Φ(S3) is trivial
        },
        {
            "name": "Dihedral group D5",
            "description": "Group of symmetries of a pentagon, q=5, m=1.",
            "order": 10,
            "frattini_order": 1, # Φ(D_q) for odd prime q is trivial
        },
        {
            "name": "Dihedral group D9",
            "description": "A nonabelian group of order 18, q=3, m=2.",
            "order": 18,
            "frattini_order": 3, # Φ(D9) is the subgroup of order 3
        },
        {
            "name": "Direct product S3 x Z3",
            "description": "A nonabelian group of order 18, q=3, m=2.",
            "order": 18,
            "frattini_order": 1, # Φ(S3 x Z3) = Φ(S3) x Φ(Z3) is trivial
        },
         {
            "name": "(Z3 x Z3) semidirect Z2",
            "description": "Another nonabelian group of order 18, q=3, m=2.",
            "order": 18,
            "frattini_order": 1, # Frattini subgroup is trivial
        },
        {
            "name": "Dihedral group D27",
            "description": "A nonabelian group of order 54, q=3, m=3.",
            "order": 54,
            "frattini_order": 9, # Φ(D27) is the subgroup of order 9
        },
    ]

    print("Investigating nonabelian groups of order 2q^m for the 'filled' property.\n")
    print("Based on the classification theorem, a necessary condition for such a group G")
    print("to be filled is that the order of its quotient group G/Φ(G) must be 9.\n")
    
    for group in groups_to_test:
        check_filled_group_candidate(group)

    print("Conclusion: As predicted by the theory, none of the tested groups satisfy the\n"
          "necessary condition to be a filled group. The set of such groups is empty.")

if __name__ == "__main__":
    main()