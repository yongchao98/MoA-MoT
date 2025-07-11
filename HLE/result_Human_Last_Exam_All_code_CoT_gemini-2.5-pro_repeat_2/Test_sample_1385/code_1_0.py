def solve_enclitic_order():
    """
    Determines and prints the correct grammatical order of a given set of Old Russian enclitics.
    """
    # The given enclitics in an arbitrary order
    given_enclitics = ["бы", "же", "еси", "бо", "мя"]

    # The established correct order of these enclitics based on linguistic rules.
    # 1. Particles (бо, же)
    # 2. Conditional particle (бы)
    # 3. Pronominal clitics (мя)
    # 4. Verbal clitics (еси)
    correct_order_reference = ["бо", "же", "бы", "мя", "еси"]

    # Sort the given list based on the index in the reference list.
    # This is a way to programmatically enforce the known grammatical order.
    sorted_enclitics = sorted(given_enclitics, key=lambda x: correct_order_reference.index(x))

    # Print the final, correctly ordered list
    print("The correct order of the enclitics is:")
    print(" -> ".join(sorted_enclitics))

solve_enclitic_order()
<<<бо, же, бы, мя, еси>>>