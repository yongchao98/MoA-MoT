def find_radical_name():
    """
    This function identifies and prints the radical name for a given Kanji.
    The information is pre-determined as there is no simple algorithmic way
    to determine a radical without a comprehensive Kanji database.
    """
    kanji_character = "謄"
    
    # The radical for 謄 is 言.
    # The Japanese name for this radical on the left side is ごんべん.
    # The Hepburn transcription is "Gonben".
    radical_name_hepburn = "Gonben"
    
    print(f"The bushu (radical) in the character '{kanji_character}' is called '{radical_name_hepburn}'.")

if __name__ == "__main__":
    find_radical_name()