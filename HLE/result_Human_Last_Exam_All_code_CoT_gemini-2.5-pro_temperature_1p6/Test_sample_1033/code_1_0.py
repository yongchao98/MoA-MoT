def solve():
    """
    This function calculates the next three terms in a letter sequence.
    It assumes the rule l3 = (l1 + l2) % 26, where l1, l2, l3 are
    0-indexed letter positions. It calculates the next three triplets
    starting with 'O'.
    """
    # The next group of letters starts with 'O'.
    l1_char = 'O'
    l1 = ord(l1_char) - ord('A')

    results = []
    # We will predict the next three terms, assuming the second letter
    # goes through 'A', 'B', 'C'.
    for l2 in range(3):
        l2_char = chr(ord('A') + l2)

        # The rule to find the third letter.
        l3 = (l1 + l2) % 26
        l3_char = chr(ord('A') + l3)

        # Form the triplet
        triplet = f"{l1_char}{l2_char}{l3_char}"
        results.append(triplet)

        # Output the equation for each triplet
        print(f"For the triplet {triplet}:")
        print(f"L1({l1_char}) = {l1}, L2({l2_char}) = {l2}")
        print(f"L3 = (L1 + L2) % 26 = ({l1} + {l2}) % 26 = {l3}, which is the letter '{l3_char}'")
        print("-" * 20)

    final_answer = " ".join(results)
    print(f"The next three capital letter sequences are: {final_answer}")
    print(f"<<<{final_answer}>>>")

solve()