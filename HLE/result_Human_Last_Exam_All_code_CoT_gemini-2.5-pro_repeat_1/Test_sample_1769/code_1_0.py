def find_hypothetical_word():
    """
    This script explains and calculates the hypothetical Modern English word for "sister"
    if the Old Norse language had not influenced Old English.
    """
    
    # The current word for sister comes from Old Norse.
    current_word = "sister"
    norse_origin = "systir"
    
    # The original, native word in Old English.
    old_english_word = "sweostor"
    
    # In the transition from Old to Middle English, vowels were simplified.
    # The 'eo' became a simple 'e' and the unstressed 'o' in the second syllable
    # was reduced.
    middle_english_word = "swester"
    
    # The Great Vowel Shift (c. 1400-1700) changed the pronunciation of long vowels.
    # The long 'e' sound in "swester" (pronounced sways-ter) would have shifted
    # to a long 'ee' sound (like in "feet").
    # The spelling would likely have adapted to reflect this.
    hypothetical_modern_word = "sweester"
    
    print("This program determines the modern word for 'sister' without Norse influence.")
    print("-" * 60)
    print(f"The current word '{current_word}' is not native to English.")
    print(f"It was adopted from the Old Norse word: '{norse_origin}'.")
    print("-" * 60)
    print("\nThe native Old English word would have evolved as follows:")
    
    # The user requested to "output each number in the final equation".
    # We will show the evolution step-by-step.
    print(f"1. The original Old English word was: {old_english_word}")
    print(f"2. Through vowel simplification, it would become the Middle English word: {middle_english_word}")
    print(f"3. After the Great Vowel Shift, the final Modern English word would be: {hypothetical_modern_word}")
    
    print("-" * 60)
    print(f"Therefore, without Norse influence, the Modern English word for 'sister' would likely be '{hypothetical_modern_word}'.")

find_hypothetical_word()