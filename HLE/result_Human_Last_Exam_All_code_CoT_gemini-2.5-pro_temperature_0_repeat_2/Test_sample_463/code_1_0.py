def solve():
    """
    This function identifies and prints the indices of taxa that perform photochemical synthesis.
    The final list of indices is constructed based on biological knowledge of each taxon.
    """
    # List of indices for taxa that undergo photochemical synthesis
    photochemical_taxa_indices = [2, 3, 5, 6, 7, 9, 10, 12]

    # The "equation" is the list of numbers themselves.
    # We will print them joined by a comma as requested.
    # The map(str, ...) converts each number in the list to a string.
    # The ",".join(...) concatenates these strings with a comma in between.
    print(",".join(map(str, photochemical_taxa_indices)))

solve()