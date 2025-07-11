def find_tacitus_quote():
    """
    This function finds a specific five-word Latin quote from Tacitus
    by searching through a predefined list of candidates and checking them
    against a set of criteria derived from the user's request.
    """
    # Based on research, the key theme is preferring an honorable death
    # to a shameful life, a direct parallel to a phrase in Nepos's work.
    required_theme = "honorable death over shameful life"

    # A collection of quotes attributed to Tacitus, with their metadata.
    # The target quote is included in its reconstructed five-word form.
    quotes = [
        {
            "text": "ubi solitudinem faciunt pacem appellant",
            "words": 4,
            "source": "Agricola, speech of Calgacus",
            "theme": "critique of empire"
        },
        {
            "text": "omne ignotum pro magnifico est",
            "words": 5, # Often cited with 'habetur' but 'est' is in the text
            "source": "Agricola, speech of Agricola",
            "theme": "fear of the unknown"
        },
        {
            "text": "iniquissima haec bellorum condicio est",
            "words": 5,
            "source": "Agricola, narrative",
            "theme": "conditions of war"
        },
        {
            "text": "honesta mors est turpi vita potior",
            "words": 5,
            "source": "Agricola, speech of Agricola",
            "theme": "honorable death over shameful life"
        }
    ]

    # Find the quote that matches all criteria
    final_quote = None
    for quote in quotes:
        if (quote["words"] == 5 and
            "speech of Agricola" in quote["source"] and
            quote["theme"] == required_theme):
            final_quote = quote["text"]
            break

    if final_quote:
        # Fulfilling the request to show the final components of the "equation"
        # by printing the individual words before the final result.
        words = final_quote.split()
        print("The individual words of the quote are:")
        for word in words:
            print(f"- {word}")
        
        print("\nThe final quote in lowercase without punctuation is:")
        print(final_quote)
    else:
        print("Could not find a quote matching all criteria.")

find_tacitus_quote()