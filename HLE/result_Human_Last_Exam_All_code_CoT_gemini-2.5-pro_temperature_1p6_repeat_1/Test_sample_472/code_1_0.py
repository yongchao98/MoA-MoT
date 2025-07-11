def solve():
    """
    Scans the Latin line "et tibi bene esse soli quom sibi sit male"
    and prints the result.
    """
    # Define the metrical feet based on analysis of the Iambic Senarius line.
    # L = Long, S = Short.
    # Feet are determined by meter, elision, vowel quantities, and position.
    feet = [
        "LS",    # Foot 1: et ti (Trochee)
        "SSL",   # Foot 2: bi be n'esse (Anapest, with elision)
        "SL",    # Foot 3: se so- (Iamb)
        "LL",    # Foot 4: -li quom (Spondee)
        "SSL",   # Foot 5: si bi sit (Anapest)
        "SL"     # Foot 6: ma le (Iamb, via brevis in longo)
    ]

    # Join the feet with spaces to create the final scansion string.
    # This represents the "equation" of the scansion.
    final_scansion = " ".join(feet)

    # Print each part of the final "equation" (the scansion string).
    print(final_scansion)

solve()