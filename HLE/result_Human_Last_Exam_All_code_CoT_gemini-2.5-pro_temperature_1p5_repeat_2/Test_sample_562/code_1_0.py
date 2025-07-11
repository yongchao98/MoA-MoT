def find_radical_name():
    """
    This function identifies and prints the name of the radical for the kanji 謄.
    """
    # The Japanese character in question.
    kanji = '謄'

    # The radical (bushu) for 謄 is the left-side component, 月.
    # This is a variant of the radical 肉 (niku, meaning 'flesh').
    # When it appears on the left side, it is known as Nikuzuki (にくづき)
    # to distinguish it from the 'moon' radical (tsukihen).
    radical_name_hepburn = "Nikuzuki"

    print(f"The bushu / radical in the Japanese character '{kanji}' is: {radical_name_hepburn}")

find_radical_name()