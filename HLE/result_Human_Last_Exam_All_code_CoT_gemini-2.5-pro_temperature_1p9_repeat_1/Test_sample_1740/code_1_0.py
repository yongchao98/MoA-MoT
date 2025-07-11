def count_let7_family_members():
    """
    This function calculates and prints the number of let-7 family members in humans.
    """
    # The canonical let-7 members (not including miR-98)
    core_let7_members = [
        "let-7a", "let-7b", "let-7c", "let-7d",
        "let-7e", "let-7f", "let-7g", "let-7i"
    ]
    
    # Other members of the let-7 family based on sequence homology
    related_members = ["miR-98"]
    
    num_core_members = len(core_let7_members)
    num_related_members = len(related_members)
    total_members = num_core_members + num_related_members
    
    print("The human let-7 family consists of several members, typically grouped into:")
    print(f"1. Core 'let-7' RNAs: {', '.join(core_let7_members)}")
    print(f"2. Related RNAs like: {', '.join(related_members)}")
    print("\nTo calculate the total:")
    print(f"Number of core members: {num_core_members}")
    print(f"Number of related members: {num_related_members}")
    print(f"Total Members = {num_core_members} + {num_related_members} = {total_members}")
    print(f"\nTo date, {total_members} members of the let-7 family have been identified in humans.")

count_let7_family_members()