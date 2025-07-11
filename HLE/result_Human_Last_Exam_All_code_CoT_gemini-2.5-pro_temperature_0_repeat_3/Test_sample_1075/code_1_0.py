def solve_legal_scenario():
    """
    This function analyzes the legal scenario involving Lewis and Marcel
    to determine the correct advice a lawyer would provide.

    The key legal points are:
    1.  The contract is for the sale of a 'good' (the painting), not a 'service'.
        Therefore, the Sale of Goods Act (SGA) applies.
    2.  Lewis made the specific purpose for the painting known to Marcel (a "centrepiece"
        for his living room), and he relied on Marcel's skill as a "highly regarded artist".
    3.  Under the SGA, this creates an implied condition that the good must be "fit for purpose".
    4.  The painting delivered (small, hastily done, of a local creek) was not what was
        described ("large", "detailed", "of Algonquin Park or Hudson Bay") and was not
        fit for its intended purpose as a centrepiece.
    5.  A breach of an implied condition allows the buyer (Lewis) to reject the goods
        and demand a full refund of the purchase price, which was $5,000.

    Answer D correctly identifies that the SGA applies and the implied condition of
    fitness for purpose was breached, entitling Lewis to a refund.
    """
    # The price of the painting
    price = 5000

    # The correct multiple-choice option
    correct_answer = 'D'

    print(f"Lewis paid Marcel ${price} for the painting.")
    print("The lawyer would advise Lewis based on the Sale of Goods Act (SGA).")
    print("The painting delivered breached the implied conditions of the SGA, including fitness for purpose.")
    print(f"Therefore, the correct option is {correct_answer}.")

solve_legal_scenario()