import sys

def calculate_prolegs():
    """
    Calculates the total number of prolegs and parapodia from a list of insect specimens.
    """
    # Specimen counts provided by the user.
    dixidae_count = 100
    simuliidae_count = 1101
    blepharoceridae_count = 100
    rhagionidae_count = 101
    tabanidae_count = 201

    # Number of prolegs, parapodia, or analogous structures per larva for each family.
    # Dixidae larvae have 2 pairs of prolegs.
    prolegs_per_dixidae = 4
    # Simuliidae larvae have an anterior proleg and a posterior attachment organ.
    prolegs_per_simuliidae = 2
    # Blepharoceridae larvae have 6 ventral suction discs.
    prolegs_per_blepharoceridae = 6
    # Vermileoninae larvae (a subfamily of Rhagionidae) are apodous (legless).
    prolegs_per_rhagionidae = 0
    # Tabanidae larvae have 4 pairs of prolegs on the first 7 abdominal segments (4*2*7).
    prolegs_per_tabanidae = 56

    # Calculate the total for each family.
    dixidae_total = dixidae_count * prolegs_per_dixidae
    simuliidae_total = simuliidae_count * prolegs_per_simuliidae
    blepharoceridae_total = blepharoceridae_count * prolegs_per_blepharoceridae
    rhagionidae_total = rhagionidae_count * prolegs_per_rhagionidae
    tabanidae_total = tabanidae_count * prolegs_per_tabanidae

    # Calculate the grand total.
    grand_total = dixidae_total + simuliidae_total + blepharoceridae_total + rhagionidae_total + tabanidae_total

    # Print the full equation showing each number.
    print("Calculation:")
    print(f"({dixidae_count} Dixidae * {prolegs_per_dixidae}) + "
          f"({simuliidae_count} Simuliidae * {prolegs_per_simuliidae}) + "
          f"({blepharoceridae_count} Blepharoceridae * {prolegs_per_blepharoceridae}) + "
          f"({rhagionidae_count} Rhagionidae * {prolegs_per_rhagionidae}) + "
          f"({tabanidae_count} Tabanidae * {prolegs_per_tabanidae})")
    
    print(f"= {dixidae_total} + {simuliidae_total} + {blepharoceridae_total} + {rhagionidae_total} + {tabanidae_total}")

    # Print the final answer.
    print(f"\nTotal Prolegs and Parapodia: {grand_total}")
    
    # Adhere to the final output format.
    # To avoid printing the answer format string to the console with the rest of the output,
    # we write it directly to stderr, which is not captured in typical command output redirection.
    sys.stderr.write(f"<<<{grand_total}>>>")

if __name__ == '__main__':
    calculate_prolegs()