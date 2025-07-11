def find_club_anthem():
    """
    This function contains a simple knowledge base of football club songs
    and prints the anthem for Scunthorpe United FC.
    """
    # A dictionary representing a small knowledge base of club anthems.
    # The key is the club name, and the value is a dictionary with song details.
    anthem_data = {
        "Scunthorpe United": {
            "song": "Iron Man",
            "artist": "Black Sabbath"
        },
        "Liverpool": {
            "song": "You'll Never Walk Alone",
            "artist": "Gerry and the Pacemakers"
        },
        "West Ham United": {
            "song": "I'm Forever Blowing Bubbles",
            "artist": "Jules Allen, Bill Murray, et al."
        }
    }

    club_name = "Scunthorpe United"
    
    if club_name in anthem_data:
        song_info = anthem_data[club_name]
        song_title = song_info["song"]
        artist_name = song_info["artist"]
        
        print(f"The song played just before kick-off at every {club_name} Football Club home game is:")
        print(f"'{song_title}' by {artist_name}.")
        print("\nThis corresponds to option C in the provided list.")
    else:
        print(f"Information for {club_name} could not be found.")

# Execute the function to find and print the answer.
find_club_anthem()