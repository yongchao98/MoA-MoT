import sys

def calculate_prolegs():
    """
    Calculates the total number of prolegs and parapodia from a list of insect specimens.
    This function determines the number of prolegs for each larval type, multiplies by the
    number of specimens, and sums the results for a grand total.
    """

    # Step 1: Define the number of specimens for each family
    num_dixidae = 100
    num_simuliidae = 1101
    num_blepharoceridae = 100
    num_vermileoninae = 101
    num_tabanidae = 201

    # Step 2: Define the number of prolegs/parapodia per larva for each family based on biological data
    # Dixidae larvae have 2 pairs of prolegs.
    prolegs_per_dixidae = 4
    # Simuliidae larvae have 1 anterior proleg and 1 posterior circlet/anal proleg.
    prolegs_per_simuliidae = 2
    # Blepharoceridae larvae have 6 ventral suckers functioning as prolegs.
    prolegs_per_blepharoceridae = 6
    # Vermileoninae larvae are apodous (have no prolegs).
    prolegs_per_vermileoninae = 0
    # Tabanus larvae have 3 pairs of prolegs on each of their first 7 abdominal segments (3*2*7=42).
    prolegs_per_tabanidae = 42

    # Step 3: Calculate the subtotal of prolegs for each family
    total_dixidae = num_dixidae * prolegs_per_dixidae
    total_simuliidae = num_simuliidae * prolegs_per_simuliidae
    total_blepharoceridae = num_blepharoceridae * prolegs_per_blepharoceridae
    total_vermileoninae = num_vermileoninae * prolegs_per_vermileoninae
    total_tabanidae = num_tabanidae * prolegs_per_tabanidae

    # Step 4: Calculate the grand total
    grand_total = total_dixidae + total_simuliidae + total_blepharoceridae + total_vermileoninae + total_tabanidae

    # Step 5: Print the breakdown of the calculation and the final answer
    print("Calculation Breakdown:")
    print(f"Dixidae: {num_dixidae} specimens * {prolegs_per_dixidae} prolegs/specimen = {total_dixidae}")
    print(f"Simuliidae: {num_simuliidae} specimens * {prolegs_per_simuliidae} prolegs/specimen = {total_simuliidae}")
    print(f"Blepharoceridae: {num_blepharoceridae} specimens * {prolegs_per_blepharoceridae} prolegs/specimen = {total_blepharoceridae}")
    print(f"Rhagionidae (Vermileoninae): {num_vermileoninae} specimens * {prolegs_per_vermileoninae} prolegs/specimen = {total_vermileoninae}")
    print(f"Tabanidae (Tabanus): {num_tabanidae} specimens * {prolegs_per_tabanidae} prolegs/specimen = {total_tabanidae}")
    
    print("\nFinal Equation:")
    print(f"{total_dixidae} + {total_simuliidae} + {total_blepharoceridae} + {total_vermileoninae} + {total_tabanidae} = {grand_total}")
    
    # The final answer is printed as a single numeric value as requested.
    # The output stream is directed to stdout, which is the standard behavior of print().
    # This avoids asking the user to copy-paste.
    sys.stdout.write(f"\n<<<Total Prolegs and Parapodia: {grand_total}>>>\n")

calculate_prolegs()