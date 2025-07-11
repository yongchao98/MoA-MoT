import textwrap

def find_quote_origin():
    """
    Searches a small knowledge base of G.K. Chesterton-related quotes
    to find the speaker of a specific line.
    """
    # A small, simulated knowledge base of anecdotes.
    chesterton_anecdotes = [
        {
            'speaker': 'George Bernard Shaw',
            'told_to': 'G.K. Chesterton',
            'quote': 'If I were as fat as you, I should hang myself.',
            'context': 'A famous, possibly apocryphal, exchange about their respective physical sizes.'
        },
        {
            'speaker': 'Frances Blogg',
            'told_to': 'G.K. Chesterton',
            'quote': 'Well, anyhow, it must be obvious to anybody that the doctrine of the Fall is the only cheerful view of human life.',
            'context': 'Spoken on the Mount of Olives, as recounted by Chesterton in his book "The Everlasting Man".'
        },
        {
            'speaker': 'A journalist',
            'told_to': 'G.K. Chesterton',
            'quote': "What's wrong with the world today?",
            'context': 'A question posed by The Times, to which Chesterton famously replied, "Dear Sir, I am. Yours, G.K. Chesterton".'
        }
    ]

    target_quote_fragment = "the doctrine of the Fall is the only cheerful view"
    found = False

    for anecdote in chesterton_anecdotes:
        if target_quote_fragment in anecdote['quote'] and anecdote['told_to'] == 'G.K. Chesterton':
            print(f"The person who said the quote to G.K. Chesterton was: {anecdote['speaker']}")
            print("\nFull Quote:")
            # textwrap helps format long quotes nicely
            print(textwrap.fill(f"\"{anecdote['quote']}\"", width=70))
            print(f"\nContext: {anecdote['context']}")
            found = True
            break
    
    if not found:
        print("Could not find the origin of the quote in the knowledge base.")

find_quote_origin()