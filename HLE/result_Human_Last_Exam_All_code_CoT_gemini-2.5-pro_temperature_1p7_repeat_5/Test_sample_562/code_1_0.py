def find_kanji_radical():
    """
    This function identifies the radical for the character 謄 and prints its name.
    """
    # The target character.
    character = '謄'
    
    # The radical for 謄 is 言.
    # When the radical 言 appears on the left side of a character, it is called "gonben" (ごんべん).
    # "Gonben" is the Hepburn transcription.
    radical_name_hepburn = "Gonben"

    print(f"The bushu (radical) in the Japanese character '{character}' is '{radical_name_hepburn}'.")

if __name__ == '__main__':
    find_kanji_radical()