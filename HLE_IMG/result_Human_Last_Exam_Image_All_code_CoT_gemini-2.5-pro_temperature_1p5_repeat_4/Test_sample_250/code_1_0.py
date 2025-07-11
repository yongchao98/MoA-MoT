def identify_music_piece():
    """
    This function identifies a musical piece based on its known characteristics.
    The provided sheet music is in the key of F minor.
    The chord progression in the first four measures is:
    F- | E-7b5 | Eb-6 | Dbmaj7 G7 |
    This progression, especially the descending bass line (F, E, Eb, Db),
    is the defining feature of the jazz standard "Yesterdays".
    """
    title = "Yesterdays"
    composer = "Jerome Kern"

    print("Identification Results:")
    print("-----------------------")
    print(f"Title: {title}")
    print(f"Composer: {composer}")
    print("-----------------------")

identify_music_piece()