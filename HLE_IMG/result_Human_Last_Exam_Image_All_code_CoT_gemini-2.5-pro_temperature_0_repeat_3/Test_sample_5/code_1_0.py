def guess_the_music():
    """
    This function identifies the piece of music from the provided sheet music
    and prints the details.
    """
    composer = "Frédéric Chopin"
    piece_title = "Piano Sonata No. 3 in B minor, Op. 58"
    movement = "IV. Finale. Presto non tanto"
    
    sonata_number = 3
    opus_number = 58

    print(f"The sheet music is from a famous piece by {composer}.")
    print(f"Title: {piece_title}")
    print(f"Movement: {movement}")
    
    print("\nAs requested, here are the numbers from the title:")
    print(f"Sonata Number: {sonata_number}")
    print(f"Opus Number: {opus_number}")

guess_the_music()