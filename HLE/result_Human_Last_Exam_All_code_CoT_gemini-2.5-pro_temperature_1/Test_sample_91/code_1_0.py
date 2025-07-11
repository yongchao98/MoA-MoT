def analyze_garner_melody():
    """
    Analyzes the musical scale used by Erroll Garner in a specific passage
    of "All My Loves Are You" and prints the findings.
    """
    song_title = "All My Loves Are You"
    artist = "Erroll Garner"
    album = "Afternoon Of An Elf (1986)"
    time_segment = "39-43 seconds"
    
    print(f"Analysis of the right-hand melody in '{song_title}' by {artist}")
    print(f"Time Segment: {time_segment}")
    print("-" * 50)
    
    key_of_song = "Ab Major"
    underlying_chord = "C7 (the V of F minor)"
    
    print(f"1. Tonal Context: The song is in {key_of_song}.")
    print(f"   The harmony in this section is based on a {underlying_chord} chord.")
    
    melody_notes = ['G', 'F', 'E', 'D', 'C', 'Bb']
    print(f"\n2. Melodic Transcription: The notes in the fast descending run are: {', '.join(melody_notes)}.")
    
    root_note = "C"
    scale_name = "C Mixolydian"
    scale_notes = ['C', 'D', 'E', 'F', 'G', 'A', 'Bb']
    
    print(f"\n3. Scale Identification: When played over a {underlying_chord} chord, these notes belong to the {scale_name} scale.")
    print(f"   The notes of the {scale_name} scale are: {', '.join(scale_notes)}.")
    print(f"   Garner's run is a descending fragment of this scale.")
    
    print("\n" + "="*50)
    final_answer = "Mixolydian Scale"
    print(f"Conclusion: The type of scale Garner plays is the {final_answer}.")
    print("="*50)

analyze_garner_melody()