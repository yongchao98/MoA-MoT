def find_imprint_differences():
    """
    Identifies and prints the different numbers on Teva 10mg and 20mg fluoxetine capsules.
    """
    # Numbers found on the 10mg capsule imprint (Teva 3105)
    imprint_10mg = {3, 1, 0, 5}

    # Numbers found on the 20mg capsule imprint (Teva 3106)
    imprint_20mg = {3, 1, 0, 6}

    # Find the numbers that are not common to both imprints
    # The symmetric_difference method returns a set with elements in either set, but not both
    different_numbers = sorted(list(imprint_10mg.symmetric_difference(imprint_20mg)))

    # --- Final Equation Output ---
    # As requested, outputting each number involved in the comparison.
    print("Numbers on the 10mg capsule:", ' '.join(map(str, sorted(list(imprint_10mg)))))
    print("Numbers on the 20mg capsule:", ' '.join(map(str, sorted(list(imprint_20mg)))))
    print("-" * 30)
    
    # Print the final result
    print("The different numbers are:", ' and '.join(map(str, different_numbers)))

find_imprint_differences()
<<<5 and 6>>>