def find_correct_triplets():
    """
    Identifies and prints the three correct triplets from the provided list.
    Each correct triplet contains a biological mode and two species that
    are directly related examples of that mode.
    """

    # Storing the correct triplets identified through biological analysis.
    # The triplets are ordered by their original numbering.
    correct_triplets = [
        "3) Batesian, Eristalini - color, Melipotini - color",
        "5) Müllerian, Heliconiini - color, Melipotini - color",
        "8) Aposematism, Danaus plexipus - wing, Cycnia tenera - tymbal"
    ]

    # Printing the identified correct triplets as the solution.
    print("The three correct triplets are:")
    for triplet in correct_triplets:
        print(triplet)

if __name__ == '__main__':
    find_correct_triplets()