def calculate_total_prolegs():
    """
    Calculates the total number of prolegs and parapodia from a given list of insect specimens.
    """
    # Step 1: Define specimen counts and proleg counts per specimen for each family.
    specimen_counts = {
        "Dixidae spp.": 100,
        "Simuliidae spp.": 1101,
        "Blepharoceridae spp.": 100,
        "Rhagionidae (Vermileoninae spp.)": 101,
        "Tabanidae (Tabanus spp.)": 201
    }

    prolegs_per_larva = {
        "Dixidae spp.": 4,
        "Simuliidae spp.": 2,
        "Blepharoceridae spp.": 6,
        "Rhagionidae (Vermileoninae spp.)": 14,
        "Tabanidae (Tabanus spp.)": 42
    }

    total_prolegs = 0
    subtotals = []
    
    print("Calculating subtotals for each family:")
    
    # Step 2: Calculate the subtotal for each family and print the breakdown.
    for family, count in specimen_counts.items():
        prolegs = prolegs_per_larva[family]
        subtotal = count * prolegs
        total_prolegs += subtotal
        subtotals.append(subtotal)
        print(f"- {family}: {count} specimens * {prolegs} prolegs/specimen = {subtotal} prolegs")

    # Step 3: Print the final summation equation and the total.
    equation_str = " + ".join(map(str, subtotals))
    print(f"\nFinal Equation: {equation_str} = {total_prolegs}")
    print("\nAnswer:")
    print(total_prolegs)

calculate_total_prolegs()