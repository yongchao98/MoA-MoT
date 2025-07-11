import json

def get_song_info():
    """
    Provides information about the song based on the provided sheet music.
    """
    song_details = {
        "Title": "Joya",
        "Musical Style": "Hard Bop / Jazz",
        "Country of Origin": "USA"
    }

    print("Title: " + song_details["Title"])
    print("Musical Style: " + song_details["Musical Style"])
    print("Country of Origin: " + song_details["Country of Origin"])

    # This part is just to fulfill the "final answer format" requirement,
    # the main answer is printed above.
    final_answer_data = {
        "title": song_details["Title"],
        "style": song_details["Musical Style"],
        "origin": song_details["Country of Origin"]
    }
    return json.dumps(final_answer_data)

# Execute the function to print the song information
song_info_json = get_song_info()
# The final answer is constructed from the data and formatted as requested.
# final_answer_string = f"<<<Title: {json.loads(song_info_json)['title']}, Musical Style: {json.loads(song_info_json)['style']}, Country of Origin: {json.loads(song_info_json)['origin']}>>>"
# print(final_answer_string)