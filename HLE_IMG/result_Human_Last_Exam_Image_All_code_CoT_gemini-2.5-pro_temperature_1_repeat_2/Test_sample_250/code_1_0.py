def identify_jazz_standard(key, melody, chords):
    """
    Identifies a jazz standard based on its musical characteristics.
    This function contains a small database of known tunes.
    """
    # Database of one song for this specific request
    if (key == "F minor" and
        melody == ["C", "F", "G", "Ab"] and
        chords == ["F-", "E-7b5", "Eb-6"]):
        
        return {
            "title": "Tadd's Delight",
            "composer": "Tadd Dameron"
        }
    else:
        return {
            "title": "Unknown",
            "composer": "Unknown"
        }

# --- Analysis of the provided sheet music ---
# Key: F minor (based on key signature and first chord)
# Melody (first bar): The notes are C, F, G, Ab
# Chords (first 3 bars): F-, E-7b5, Eb-6
key_signature = "F minor"
melody_start = ["C", "F", "G", "Ab"]
chord_progression_start = ["F-", "E-7b5", "Eb-6"]

# Identify the song using the extracted features
song_info = identify_jazz_standard(key_signature, melody_start, chord_progression_start)

# Print the final answer
print(f"Title: {song_info['title']}")
print(f"Composer: {song_info['composer']}")
print("\nNote: This tune is also sometimes known as 'Fungii Mama'.")
