def solve_mimicry_puzzle():
    """
    This function identifies and prints the three correct triplets from the provided list,
    ordered by their original numbering.
    """
    # Based on biological analysis, the following triplets are identified as correct.
    correct_triplets = {
        4: "Gilbertian, Passiflora - leaf, Myrmecia chrysogaster - venom",
        5: "Müllerian, Heliconiini - color, Melipotini - color",
        8: "Aposematism, Danaus plexipus - wing, Cycnia tenera - tymbal"
    }

    # Print the triplets in ascending order by their number.
    print("The three correct triplets in ascending order are:")
    for number in sorted(correct_triplets.keys()):
        print(f"{number}) {correct_triplets[number]}")

solve_mimicry_puzzle()