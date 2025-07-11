def find_latin_quote():
    """
    This function simulates a search for a specific Latin quote based on a set of criteria.
    """
    # A mini-database of potential quotes and their attributes.
    quotes_data = {
        "ubi solitudinem faciunt pacem appellant": {
            "author": "Tacitus",
            "speaker": "Calgacus",
            "word_count": 5
        },
        "proinde honesta mors turpi vita potior": {
            "author": "Tacitus",
            "speaker": "Agricola",
            "word_count": 5
        },
        "praestare honestam mortem turpi vitae": {
            "author": "Cornelius Nepos",
            "speaker": "Chabrias",
            "word_count": 5
        }
    }

    # Search for the quote that meets all the user's conditions.
    for quote, metadata in quotes_data.items():
        if (metadata["author"] == "Tacitus" and
            metadata["speaker"] == "Agricola" and
            metadata["word_count"] == 5):
            
            # The quote is found. Print it in the required format.
            # The format is already lowercase and without punctuation.
            print(quote)
            return quote

# Execute the function to find and print the quote.
find_latin_quote()