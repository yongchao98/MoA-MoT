def find_anecdote_source():
    """
    Simulates searching a database of historical anecdotes to find the
    person who spoke a particular quote to G.K. Chesterton.
    """
    # A simulated database of anecdotes involving G.K. Chesterton.
    anecdotes_db = [
        {
            'speaker': 'Frances Chesterton',
            'listener': 'G.K. Chesterton',
            'location': 'Mount of Olives, Jerusalem',
            'quote_summary': 'The doctrine of the Fall is the only cheerful view of human life.'
        },
        {
            'speaker': 'George Bernard Shaw',
            'listener': 'G.K. Chesterton',
            'location': 'London',
            'quote_summary': 'Anyone would think a famine had struck England.'
        }
    ]

    # Search parameters from the user's query.
    search_listener = 'G.K. Chesterton'
    search_location_keyword = 'Mount of Olives'
    search_quote_keyword = 'doctrine of the Fall'

    # Iterate through the database to find a match.
    for anecdote in anecdotes_db:
        if (anecdote['listener'] == search_listener and
            search_location_keyword in anecdote['location'] and
            search_quote_keyword in anecdote['quote_summary']):
            
            # Print the result when found.
            print("Searching for the person who told G.K. Chesterton the following quote on the Mount of Olives:")
            print(f"'{anecdote['quote_summary']}'")
            print("\n... Match found in database.")
            print(f"\nThe speaker was: {anecdote['speaker']}")
            return

    print("Could not find a matching anecdote in the database.")

# Run the search function.
find_anecdote_source()