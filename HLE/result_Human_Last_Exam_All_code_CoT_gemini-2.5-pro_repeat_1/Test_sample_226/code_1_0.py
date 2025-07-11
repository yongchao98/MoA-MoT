def calculate_total_prolegs():
    """
    Calculates the total number of prolegs and parapodia from a collection of insect specimens.
    """
    # Specimen counts for each family
    dixidae_count = 100
    simuliidae_count = 1101
    belpharoceridae_count = 100
    rhagionidae_count = 101
    tabanidae_count = 201

    # Number of prolegs/parapodia per specimen for each family
    # Dixidae larvae have 3 pairs of prolegs (6 total)
    dixidae_prolegs = 6
    # Simuliidae larvae have 1 anterior and 1 posterior proleg (2 total)
    simuliidae_prolegs = 2
    # Belpharoceridae larvae have suckers, not prolegs (0 total)
    belpharoceridae_prolegs = 0
    # Rhagionidae (Vermileoninae) larvae are legless (0 total)
    rhagionidae_prolegs = 0
    # Tabanidae larvae have rings of prolegs on 7 segments (8 prolegs * 7 segments = 56 total)
    tabanidae_prolegs = 56

    # Calculate sub-totals for each family
    total_dixidae = dixidae_count * dixidae_prolegs
    total_simuliidae = simuliidae_count * simuliidae_prolegs
    total_belpharoceridae = belpharoceridae_count * belpharoceridae_prolegs
    total_rhagionidae = rhagionidae_count * rhagionidae_prolegs
    total_tabanidae = tabanidae_count * tabanidae_prolegs

    # Calculate the grand total
    grand_total = total_dixidae + total_simuliidae + total_belpharoceridae + total_rhagionidae + total_tabanidae
    
    # Print the explanation and the final equation
    print("Calculating the total number of prolegs and parapodia:")
    print(f"Dixidae: {dixidae_count} specimens * {dixidae_prolegs} prolegs/specimen = {total_dixidae}")
    print(f"Simuliidae: {simuliidae_count} specimens * {simuliidae_prolegs} prolegs/specimen = {total_simuliidae}")
    print(f"Belpharoceridae: {belpharoceridae_count} specimens * {belpharoceridae_prolegs} prolegs/specimen = {total_belpharoceridae}")
    print(f"Rhagionidae: {rhagionidae_count} specimens * {rhagionidae_prolegs} prolegs/specimen = {total_rhagionidae}")
    print(f"Tabanidae: {tabanidae_count} specimens * {tabanidae_prolegs} prolegs/specimen = {total_tabanidae}")
    
    print("\nFinal Equation:")
    print(f"({dixidae_count} * {dixidae_prolegs}) + ({simuliidae_count} * {simuliidae_prolegs}) + ({belpharoceridae_count} * {belpharoceridae_prolegs}) + ({rhagionidae_count} * {rhagionidae_prolegs}) + ({tabanidae_count} * {tabanidae_prolegs}) = {grand_total}")

    print(f"\nTotal prolegs and parapodia in the dish: {grand_total}")

if __name__ == "__main__":
    calculate_total_prolegs()