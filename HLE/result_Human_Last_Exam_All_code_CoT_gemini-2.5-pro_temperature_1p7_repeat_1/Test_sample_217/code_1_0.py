def solve_mimicry_puzzle():
    """
    Identifies and prints the three correct triplets from the provided list.
    """
    correct_triplets = [
        "3) Batesian, Eristalini - color, Melipotini - color",
        "5) Müllerian, Heliconiini - color, Melipotini - color",
        "8) Aposematism, Danaus plexipus - wing, Cycnia tenera - tymbal"
    ]

    print("The three correct triplets are:")
    for triplet in correct_triplets:
        print(triplet)

solve_mimicry_puzzle()