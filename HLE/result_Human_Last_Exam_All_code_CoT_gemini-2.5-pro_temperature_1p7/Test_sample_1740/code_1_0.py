def count_let7_family_members():
    """
    This function counts and displays the number of canonical members
    of the let-7 microRNA family identified in humans.
    """
    # The canonical members of the human let-7 family
    let7_family_members = [
        "let-7a", "let-7b", "let-7c", "let-7d", "let-7e",
        "let-7f", "let-7g", "let-7i", "miR-98"
    ]

    count = len(let7_family_members)

    # Build and print the equation string
    equation_parts = []
    for member in let7_family_members:
        equation_parts.append(f"1")

    equation_str = " + ".join(equation_parts)

    print("The number of human let-7 family members can be calculated as:")
    print(f"{equation_str} = {count}")
    print("\nThese members are:")
    for member in let7_family_members:
        print(f"- {member}")
    
    print(f"\nTo date, {count} members of the let-7 family have been identified in humans.")

count_let7_family_members()