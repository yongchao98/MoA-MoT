def check_persistence(Ne, s):
    """
    Checks if a genomic fragment is likely to persist based on
    effective population size (Ne) and selection coefficient (s).
    Selection is considered efficient if s > 1/Ne.
    """
    # The strength of genetic drift is inversely proportional to population size.
    drift_strength = 1 / Ne

    # Determine if selection is strong enough to overcome drift.
    is_efficient = s > drift_strength
    
    print(f"--- Scenario: Effective Population Size (Ne) = {Ne} ---")
    print(f"The fitness advantage of the fragment (s) is: {s}")
    print(f"The approximate strength of genetic drift (1/Ne) is: {drift_strength:.6f}")
    
    # Print the "equation" or comparison
    print(f"Is s > 1/Ne?  {s} > {drift_strength:.6f} ?  {is_efficient}")

    if is_efficient:
        print("Result: Selection is efficient. The fragment is likely to PERSIST.\n")
    else:
        print("Result: Selection is inefficient; drift dominates. The fragment is likely to be LOST.\n")

# A small, but positive, fitness advantage for keeping the fragment.
selection_advantage = 0.0005

# Scenario 1: A small population of endosymbionts (e.g., Ne = 5,000)
# Here, selection is often inefficient.
ne_small = 5000
check_persistence(ne_small, selection_advantage)

# Scenario 2: A large population of free-living bacteria (e.g., Ne = 5,000,000)
# Here, selection is much more efficient.
ne_large = 5000000
check_persistence(ne_large, selection_advantage)