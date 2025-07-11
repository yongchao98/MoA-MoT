def find_worst_suited_jack():
    """
    Identifies the worst suited jack to open from the button in a standard 100BB game.

    In poker, the value of a suited hand is determined by its high card value,
    flush potential, and straight potential (connectivity). For suited Jacks (JXs),
    the Jack and the flush potential are constant. The value therefore decreases
    as the rank of the second card ('kicker') goes down, because both the kicker
    value and the hand's ability to make a straight are reduced.

    J2s has the lowest kicker and the worst connectivity, making it the weakest,
    or "worst," suited Jack that is still a standard open-raise from the button.
    """

    # List of possible kickers for a suited Jack, in order from best to worst
    # (excluding Ace, King, Queen, which are premium hands).
    # Ten down to Two.
    kickers = ['T', '9', '8', '7', '6', '5', '4', '3', '2']

    # The worst kicker is the last one in the list.
    worst_kicker = kickers[-1]

    # The worst suited jack to open is J with the worst kicker.
    worst_hand = f"J{worst_kicker}s"

    print(f"The suited Jack hands (JXs) are valued based on their second card.")
    print(f"A lower second card means lower kicker value and worse straight potential.")
    print(f"The hand with the lowest value that is still a standard open from the button is J2s.")
    print(f"Therefore, the worst suited jack you should open from the button is: {worst_hand}")

find_worst_suited_jack()