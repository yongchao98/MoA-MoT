import sys

def find_alt_sister():
    """
    This function explains the linguistic path from Old English 'sweostor'
    to its hypothetical modern form without Norse influence.
    """
    # The current word "sister" comes from the influence of the Old Norse word "systir".
    # To find what the word would have been, we must start with the original Old English word.
    old_english_word = "sweostor"

    print(f"The original Old English (Anglo-Saxon) word for 'sister' was: {old_english_word}")
    print("-" * 50)
    print("Without the influence of Old Norse 'systir', this word would have undergone native English sound changes.")
    print("\nHere is the likely step-by-step evolution:\n")

    # Step 1: The initial 'sw-' sound cluster often simplified to just 's-'.
    # For example, Old English 'swa' became 'so'.
    step1_word = "sostor"
    print(f"Step 1: The 'sw' sound in '{old_english_word}' simplifies to 's'.")
    print(f"   - This is similar to how Old English 'swa' became modern 'so'.")
    print(f"   - Hypothetical form: {step1_word}\n")


    # Step 2: The Old English vowel 'eo' often changed to 'u' in Middle English.
    final_word = "suster"
    print(f"Step 2: The 'eo' vowel sound shifts to 'u'.")
    print(f"   - This was a common sound change from Old to Middle English.")
    print(f"   - Applying this to '{step1_word}' gives the final form.\n")

    # Final Result and comparison with related languages.
    print("-" * 50)
    print("The final equation of this linguistic change is:")
    print(f"'{old_english_word}' (Old English) -> '{final_word}' (Hypothetical Modern English)\n")

    print("This conclusion is supported by looking at English's closest relatives:")
    print("- Dutch: zuster")
    print("- West Frisian: suster")
    print("\nTherefore, without Norse influence, the Modern English word would most likely be:")
    print(f"'{final_word}'")


if __name__ == '__main__':
    find_alt_sister()
    # The final answer is wrapped below for clarity.
    # The name of the file is not relevant, the output is what matters.
    # To run this code, you would save it as a .py file (e.g., linguistic_analysis.py)
    # and execute it from your terminal with: python linguistic_analysis.py
    sys.stdout.write("\n<<<suster>>>\n")