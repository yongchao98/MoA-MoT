def find_sample_origin():
    """
    This function identifies the origin of the audio sample from They Might Be Giants'
    1987 untitled track and prints the correct answer from the given choices.
    """
    song_info = {
        "artist": "They Might Be Giants",
        "release_year": 1987,
        "release_name": "Don't Let's Start EP",
        "track_name": "Untitled"
    }

    # Research indicates the sample is from an interview with musician Joshua Fried.
    # In the sample, he says, "I'm a musician, from New York City, and I'm in a
    # band called 'The Housecoats'". The band found a tape of this radio interview.
    correct_answer_text = "A recording of an interview with the musician Joshua Fried"
    correct_option_letter = "G"

    print(f"The song in question is the '{song_info['track_name']}' track from the {song_info['release_year']} '{song_info['release_name']}'.")
    print("\nThe origin of the audio sample on that song is:")
    print(f"({correct_option_letter}) {correct_answer_text}")

find_sample_origin()