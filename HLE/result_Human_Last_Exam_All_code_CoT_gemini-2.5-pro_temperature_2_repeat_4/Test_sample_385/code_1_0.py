def find_latin_quote():
    """
    This function finds and formats the requested five-word Latin quote.
    The quote from Tacitus's "Agricola" (chapter 34) is "Praestat honesta mors turpi vita".
    """
    # The five-word Latin quote attributed to Agricola
    quote = "Praestat honesta mors turpi vita"
    
    # Convert to lowercase and ensure no punctuation
    formatted_quote = quote.lower().strip()
    
    # Print the final formatted quote
    print(formatted_quote)

find_latin_quote()