def find_tacitus_quote():
    """
    This function stores and prints the specific five-word Latin quote requested.
    
    The quote is: "omnia praesumi pro hoste deserta"
    - It is five words long.
    - It is from the author Tacitus (Historiae 1.63).
    - It means "all things abandoned are presumed to be for the enemy," which
      closely resembles the strategic advice in Nepos' Life of Chabrias, Ch. 4.
    
    The output will be in all lowercase, without punctuation.
    """
    
    quote_word_1 = "omnia"
    quote_word_2 = "praesumi"
    quote_word_3 = "pro"
    quote_word_4 = "hoste"
    quote_word_5 = "deserta"
    
    final_quote = f"{quote_word_1} {quote_word_2} {quote_word_3} {quote_word_4} {quote_word_5}"
    
    print(final_quote)

find_tacitus_quote()