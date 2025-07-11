def solve_runestone_id():
    """
    Identifies the Ingvar runestone from the provided image fragment.
    """
    # Step 1 & 2: Analyze the runes in the image fragment.
    # The runes are from the Younger Futhark alphabet.

    # Top row transliteration: ...ᚦᚢᚱᛒᛁᚢᚱᚾ... which is "thurbiurn", a spelling of the name Þorbjörn (Thorbjorn).
    # Middle row transliteration: ...ᛋᚢᚾ + ᛋᛁᚾ... which is "sun sin" meaning "his son".
    # Bottom row transliteration: ...ᚼᛅᚾ... which is "han" meaning "he".

    fragment_transliteration = {
        "Top Row": "...thurbiurn...",
        "Middle Row": "...sun sin...",
        "Bottom Row": "...han..."
    }

    # Step 3: Match the fragment to a known Ingvar runestone inscription.
    # This specific sequence of words is found on the Sö 281 runestone.
    # The relevant part of the Sö 281 inscription is:
    # Old Norse: "... at Þorbiorn, sun sinn. Hann ændaðis ..."
    # English: "... in memory of Þorbjörn, his son. He met his end ..."

    runestone_id = "Sö 281"
    
    # Step 4 & 5: Print the identification and final answer.
    print("The runestone in the image has been identified based on its inscription.")
    print(f"Fragment Analysis:")
    for row, text in fragment_transliteration.items():
        print(f"- {row}: Translates to '{text}'")

    print(f"\nThis fragment matches the inscription on the Ingvar runestone with the ID: {runestone_id}")
    
    # As requested, printing each number in the final ID.
    id_number = "281"
    print(f"The number in the ID is {id_number}. The individual digits are:")
    for digit in id_number:
        print(digit)

solve_runestone_id()