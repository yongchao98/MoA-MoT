def find_hypothetical_word():
    """
    This function explains and calculates the hypothetical Modern English word
    for "sister" without Norse influence.
    """
    old_norse_word = "systir"
    modern_word = "sister"
    old_english_word = "sweostor"

    print(f"The modern English word '{modern_word}' is derived from the Old Norse word '{old_norse_word}'.")
    print(f"Without the Viking invasions, the native Old English word '{old_english_word}' would have evolved naturally.")
    print("\nHere is its likely linguistic path:")
    print(f"1. Old English (c. 450-1150): The word was 'sweostor'.")

    # Linguistic change 1: Simplification of the 'eo' diphthong
    # This is analogous to OE 'seofon' -> ME 'seven'
    middle_english_word = "swester"
    print(f"2. Middle English (c. 1150-1500): The 'eo' sound typically simplified to 'e'.")
    print(f"   Therefore, '{old_english_word}' would likely have become '{middle_english_word}'.")

    # Linguistic change 2: Stability into Modern English
    # The form 'swester' is phonetically stable and would likely not undergo major changes.
    hypothetical_modern_word = "swester"
    print(f"3. Modern English (c. 1500-present): The form '{middle_english_word}' would have been very stable.")
    print(f"   It would likely have carried into Modern English as '{hypothetical_modern_word}'.")
    print("\nThis word, or a close variant, still exists in some regional English and Scots dialects.")

    print("\n--------------------------------------------------")
    print(f"The hypothetical Modern English word for 'sister' would be:")
    print(hypothetical_modern_word)
    print("--------------------------------------------------")

    return hypothetical_modern_word

# Execute the function to find and print the answer.
final_answer = find_hypothetical_word()
# The final answer is returned in the specified format below.
# For the purpose of this script, we print it here as well.
# print(f"\n<<<{final_answer}>>>")
