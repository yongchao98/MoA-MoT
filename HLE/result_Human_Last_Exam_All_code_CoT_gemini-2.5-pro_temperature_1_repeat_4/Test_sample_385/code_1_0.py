def find_tacitus_quote():
    """
    This function stores the words of the specified Latin quote,
    formats them as requested, and prints the result.
    """
    # The five words of the Latin quote from Tacitus' Agricola (Chapter 33).
    quote_words = ["honesta", "mors", "turpi", "vita", "potior"]

    # Join the words with spaces and convert to lowercase.
    final_quote = " ".join(quote_words).lower()

    # Print the final formatted quote.
    print(final_quote)

find_tacitus_quote()