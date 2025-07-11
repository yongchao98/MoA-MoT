def find_hypothetical_word():
    """
    This function explains the linguistic reasoning to find what the word "sister"
    would be in Modern English without Norse influence.
    """
    
    # The current word "sister" comes from Old Norse "systir",
    # introduced during the Viking invasions.
    old_norse_word = "systir"
    
    # The original, native Old English word was "sweostor".
    old_english_word = "sweostor"
    
    # To see how "sweostor" might have evolved, we can look at other Old English
    # family terms and their modern equivalents.
    # OE brōþor -> ModE brother
    # OE mōdor  -> ModE mother
    # This shows a common pattern where the unstressed "-or" ending becomes "-er".
    
    # Applying this pattern to "sweostor" would give us "swester".
    # This form is supported by two key pieces of evidence:
    # 1. The word "swester" did exist in Middle English, alongside "sister", before "sister" won out.
    # 2. The cognate word in Modern German, which did not undergo the same Norse influence, is "Schwester".
    
    hypothetical_word = "swester"

    print(f"The current word 'sister' comes from the Old Norse word '{old_norse_word}'.")
    print(f"However, the original Old English word was '{old_english_word}'.")
    print("\nIf this word had evolved naturally into Modern English, following the same pattern as 'brother' and 'mother',")
    print("it would have become:")
    print(f"\n---> {hypothetical_word} <---")

find_hypothetical_word()