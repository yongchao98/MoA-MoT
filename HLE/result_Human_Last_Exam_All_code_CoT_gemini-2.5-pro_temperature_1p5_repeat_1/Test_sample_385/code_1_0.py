def solve_latin_quote_puzzle():
    """
    This function solves the puzzle about a five-word Latin quote from Tacitus.

    1. The quote is from Tacitus's "Agricola," in a speech before the Battle of Mons Graupius.
    2. It's attributed to Agricola and closely resembles a phrase in Cornelius Nepos's "Life of Chabrias," which discusses honorable death over shameful life.
    3. The core phrase in Tacitus is "honesta mors turpi vita potior" ("an honorable death is preferable to a shameful life"). This is four words.
    4. To meet the user's requirement for a five-word quote, the commonly omitted Latin verb "est" ("is") must be included.
    5. The final five-word phrase is then "honesta mors turpi vita potior est".
    6. The output is formatted in all lowercase and without punctuation as requested.
    """
    quote = "honesta mors turpi vita potior est"
    print(quote)

solve_latin_quote_puzzle()