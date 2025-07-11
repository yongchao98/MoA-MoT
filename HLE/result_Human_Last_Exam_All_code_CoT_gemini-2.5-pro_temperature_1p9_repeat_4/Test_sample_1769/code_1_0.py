import sys
import io

# Buffer to capture print output
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

def derive_word():
    """
    This function simulates the linguistic evolution of the Old English word for "sister"
    in a hypothetical timeline without Norse influence.
    """
    # Step 1: Start with the original Old English word.
    # The Old English word for sister was 'sweostor'. This is our starting point.
    old_english_word = "sweostor"
    print(f"The starting Old English word is: {old_english_word}")

    # Step 2: Simulate the first sound change: Simplification of the initial consonant cluster.
    # The initial 'sw-' sound often simplified to just 's-' during the Middle English period.
    # Transformation: sweostor -> seostor
    middle_english_form_1 = "s" + old_english_word[2:]
    print(f"After Step 1 (sw -> s): {middle_english_form_1}")

    # Step 3: Simulate the second sound change: Vowel shift.
    # The 'eo' diphthong in 'seostor' would likely shift to 'u'. This is based on
    # 'suster' being a very common Middle English form of the word.
    # Transformation: seostor -> sustor
    middle_english_form_2 = middle_english_form_1.replace("eo", "u")
    print(f"After Step 2 (eo -> u): {middle_english_form_2}")
    
    # Step 4: Simulate the third sound change: Suffix leveling.
    # The ending '-or' often leveled to '-er' in kinship terms, by analogy with
    # words like 'brother' and 'mother' (from OE 'broþor' and 'modor').
    # Transformation: sustor -> suster
    hypothetical_word = middle_english_form_2.replace("or", "er")
    print(f"After Step 3 (-or -> -er): {hypothetical_word}")

    print("\n---")
    print("Conclusion:")
    print(f"The Old English '{old_english_word}' would likely have become Middle English '{hypothetical_word}'.")
    print("In Modern English, this would probably be spelled 'Suster',")
    print("with the 'u' pronounced like the vowel in 'custom' or 'putt'.")
    print("\nFinal Hypothetical Word:")
    print(hypothetical_word.capitalize())

derive_word()

# Restore stdout and get the captured output
sys.stdout = old_stdout
output = captured_output.getvalue()

# Final display
print(output)