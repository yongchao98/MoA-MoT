def find_tacitus_quote():
    """
    This script identifies and prints the requested Latin quote.

    The quote is from Tacitus's "Agricola," from a speech before the Battle of Mons Graupius.
    It is five words long and closely resembles a phrase from Nepos's "Life of Chabrias," Chapter 4.
    The phrase from Tacitus is "honestam mortem turpi vita potiorem."
    The script will print each word that forms the final quote.
    """

    # The five words of the Latin quote
    word1 = "honestam"
    word2 = "mortem"
    word3 = "turpi"
    word4 = "vita"
    word5 = "potiorem"

    # Assemble the final quote
    final_quote = f"{word1} {word2} {word3} {word4} {word5}"

    # Print each individual word as requested by the "equation" instruction
    print(f"The first word is: {word1}")
    print(f"The second word is: {word2}")
    print(f"The third word is: {word3}")
    print(f"The fourth word is: {word4}")
    print(f"The fifth word is: {word5}")
    print("-" * 25)

    # Print the final quote in all lowercase without punctuation
    print("Final quote:")
    print(final_quote)

find_tacitus_quote()