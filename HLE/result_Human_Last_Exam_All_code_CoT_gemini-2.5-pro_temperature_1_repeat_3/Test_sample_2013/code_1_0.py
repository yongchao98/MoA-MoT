def find_musical_joke_measures():
    """
    This function identifies and formats the answer to the musical joke query.
    The famous joke occurs in Haydn's Piano Concerto No. 11 in D major.
    The cello section's flourish begins in measure 161.
    They are subsequently off-beat for measures 163 and 164.
    They recover in measure 165 by skipping a beat and re-entering correctly.
    The full range of the joke, from the start of the flourish to the completed recovery,
    is measures 161-165.
    """
    composer_surname = "Haydn"
    # The Hoboken-Verzeichnis (Hob.) is the standard catalog for Haydn's works.
    catalog_number = "Hob. XVIII:11"
    start_measure = 161
    end_measure = 165
    
    # Format the range as a hyphenated string
    measure_range = f"{start_measure}-{end_measure}"
    
    # Format the final answer string
    final_answer = f"{composer_surname}, {catalog_number}, {measure_range}"
    
    print(final_answer)

find_musical_joke_measures()