def find_hypothetical_word():
    """
    This function demonstrates the likely linguistic evolution of the Old English
    word for "sister" into a hypothetical Modern English equivalent, had it not
    been replaced by a word from Old Norse.
    """
    old_english_word = "sweostor"
    
    print(f"The native Old English word for 'sister' was: {old_english_word}")
    
    # Step 1: The 'eo' diphthong often simplified to 'e'.
    # This is a common sound change from Old to Middle English.
    middle_english_step1 = old_english_word.replace("eo", "e")
    print(f"Step 1: The vowel 'eo' would likely simplify to 'e', resulting in: {middle_english_step1}")

    # Step 2: The unstressed ending '-or' in kinship terms like 'broþor' and 'modor'
    # became '-er' in Modern English ('brother', 'mother').
    hypothetical_modern_english_word = middle_english_step1[:-2] + "er"
    print(f"Step 2: The ending '-or' would likely shift to '-er', resulting in: {hypothetical_modern_english_word}")
    
    print("\n---")
    print(f"Therefore, without the Norse influence that gave us 'sister', the Modern English word would likely have been '{hypothetical_modern_english_word}'.")
    print(f"This is supported by the survival of 'suster' in some dialects and the Scots language.")

find_hypothetical_word()