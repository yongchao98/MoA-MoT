def calculate_total_prolegs():
    """
    Calculates the total number of prolegs and parapodia from a collection of insect specimens.
    """
    # Specimen counts for each family
    specimen_counts = {
        "Dixidae": 100,
        "Simuliidae": 1101,
        "Belpharoceridae": 100,
        "Rhagionidae": 101,
        "Tabanidae": 201
    }

    # Number of prolegs/parapodia per larva for each family
    prolegs_per_larva = {
        "Dixidae": 4,           # 2 pairs of prolegs
        "Simuliidae": 2,        # 1 anterior and 1 posterior proleg
        "Belpharoceridae": 0,   # Ventral suckers are not prolegs/parapodia
        "Rhagionidae": 16,      # Paired welts on 8 segments (2 * 8)
        "Tabanidae": 42         # 3 pairs of prolegs on 7 segments (3 * 2 * 7)
    }

    # Calculate the subtotal of prolegs for each family
    subtotals = {family: count * prolegs_per_larva[family] for family, count in specimen_counts.items()}

    # Calculate the grand total
    grand_total = sum(subtotals.values())

    # Print the breakdown of the calculation
    print("Calculating the total number of prolegs and parapodia:")
    for family, count in specimen_counts.items():
        prolegs = prolegs_per_larva[family]
        subtotal = subtotals[family]
        print(f"{family}: {count} specimens * {prolegs} prolegs/specimen = {subtotal} prolegs")

    # Print the final equation
    equation_parts = [str(s) for s in subtotals.values()]
    equation_str = " + ".join(equation_parts)
    print(f"\nTotal Calculation:\n{equation_str} = {grand_total}")

    print(f"\nThe total number of prolegs and parapodia in the dish is: {grand_total}")

calculate_total_prolegs()