def find_tacitus_quote():
    """
    This function stores and prints the requested five-word Latin quote.

    The quote is from Tacitus' "Agricola," Chapter 34, from a speech
    Agricola gives before the Battle of Mons Graupius.
    """
    # The five-word Latin quote attributed to Agricola
    quote_words = ["maiores", "vestros", "et", "posteros", "cogitate"]

    # The final quote is formed by joining the words with spaces.
    final_quote = " ".join(quote_words)

    # Print the quote in lowercase as requested.
    print(final_quote.lower())

find_tacitus_quote()