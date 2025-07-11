def find_highest_recombinant_frequency_location():
    """
    This function models the E. coli interrupted mating experiment to find
    the location of the highest recombinant frequency.
    """

    # Step 1: Define the gene markers and their time of entry in minutes.
    # The problem states that thr+ is transferred first, followed by azy,
    # and the gene order is thr-azi-gal. This means thr is closest to the
    # origin of transfer, followed by azi, and then gal.
    # We assign arbitrary entry times to reflect this order.
    gene_entry_times = {
        'thr': 8,
        'azi': 15,
        'gal': 25
    }

    # Step 2: The principle of interrupted mating states that genes closer
    # to the origin of transfer (with lower entry times) are transferred
    # more frequently. Higher transfer frequency leads to a higher
    # frequency of recombination.
    
    # We find the gene with the minimum entry time.
    location_highest_freq = min(gene_entry_times, key=gene_entry_times.get)
    min_time = gene_entry_times[location_highest_freq]
    
    # We can also explicitly state the gene order based on sorted times.
    gene_order_sorted = sorted(gene_entry_times, key=gene_entry_times.get)

    # Step 3: Print the results and the reasoning.
    print(f"Assumed entry times (minutes): thr = {gene_entry_times['thr']}, azi = {gene_entry_times['azi']}, gal = {gene_entry_times['gal']}")
    print(f"The determined order of gene transfer is: {' -> '.join(gene_order_sorted)}")
    print(f"The gene with the lowest time of entry is '{location_highest_freq}' at {min_time} minutes.")
    print("\nConclusion:")
    print(f"The highest frequency of recombinants occurs at the locus of the first gene transferred, which is '{location_highest_freq}'.")
    print("This corresponds to the location immediately after the 'thr+' marker enters the recipient cell.")

find_highest_recombinant_frequency_location()