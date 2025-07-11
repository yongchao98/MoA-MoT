import requests
import string

def find_black_player():
    """
    Finds the black player of a chess game by querying the Lichess Masters database API.
    """
    # The PGN of the game provided by the user. We only need the first 15 moves for a unique match.
    pgn_moves = "1.d4 Nf6 2. Nf3 d5 3. g3 e6 4. Bg2 Be7 5. O-O O-O 6. b3 c5 7. dxc5 Bxc5 8. c4 dxc4 9. Qc2 Qe7 10. Nbd2 Nc6 11. Nxc4 b5 12. Nce5 Nb4 13. Qb2 Bb7 14. a3 Nc6 15. Nd3 Bb6"

    # The list of possible players.
    choices = [
        "Anand, Viswanathan", "Karpov, Anatoly", "Keymer, Vincent", "Ding, Liren",
        "Aronian, Levon", "Radjabov, Teimour", "Kasparov, Garry", "Firouzja, Alireza",
        "So, Wesley", "Giri, Anish", "Nepomniachtchi, Ian", "Kramnik Vladimir",
        "Fischer, Robert", "Grischuck Alexander", "Niemann, Hans Moke",
        "Caruana, Fabiano", "Carlsen, Magnus"
    ]
    
    # Create a mapping from letters (A, B, C...) to player names.
    choice_map = {letter: name for letter, name in zip(string.ascii_uppercase, choices)}

    # The API endpoint for the Lichess Masters database.
    api_url = "https://explorer.lichess.ovh/master"

    # The parameters for the GET request. We send the PGN moves.
    params = {'pgn': pgn_moves}

    try:
        # Make the request to the Lichess API
        response = requests.get(api_url, params=params)
        response.raise_for_status()  # Raise an exception for bad status codes (4xx or 5xx)

        # Parse the JSON response
        data = response.json()

        # Check if any games were found
        if 'topGames' in data and data['topGames']:
            # The first game in the list is the most relevant match.
            top_game = data['topGames'][0]
            
            # Extract the name of the player of the black pieces.
            black_player_name = top_game['black']['name']
            
            # Find the matching name in our choices
            for letter, name in choice_map.items():
                # Lichess may use a slightly different name format (e.g. adding middle initials)
                # so we check if the name from the API starts with one of our choices.
                # The names are in "Lastname, Firstname" format.
                if black_player_name.startswith(name):
                    print(f"The player of the black pieces was: {black_player_name}")
                    print(f"This corresponds to option {letter}: {name}")
                    print(f"<<<{letter}>>>")
                    return
            
            print(f"Found player '{black_player_name}', but could not match to the choices.")

        else:
            print("Could not find a matching game in the database.")

    except requests.exceptions.RequestException as e:
        print(f"An error occurred while contacting the Lichess API: {e}")
    except KeyError:
        print("Received an unexpected format from the API.")

if __name__ == '__main__':
    find_black_player()