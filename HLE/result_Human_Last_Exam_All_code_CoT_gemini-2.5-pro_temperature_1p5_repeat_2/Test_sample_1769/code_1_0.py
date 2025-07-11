def find_alternate_word():
    """
    Calculates and prints the hypothesized Modern English word for 'sister'
    if the Norse had never influenced the language.
    """
    old_english = "sweostor"
    middle_english = "suster"
    modern_english = "suster"

    print("This script traces the linguistic evolution of the Old English word for 'sister' to its potential modern form without Norse influence.")
    print("-" * 50)
    print(f"1. The current word 'sister' comes from Old Norse 'systir', which replaced the native Old English word: '{old_english}'.")
    print(f"2. Without that influence, '{old_english}' would have evolved naturally. It was often simplified in Middle English to '{middle_english}'.")
    print(f"3. The vowel 'u' in '{middle_english}' would likely have undergone the 'foot-strut split' sound change (like in 'sunne' -> 'sun'), leading to the modern form.")
    print("-" * 50)
    print("The final derivation path (the 'equation') is:")
    # The prompt requests that we "output each number in the final equation",
    # which is interpreted here as showing each word in the evolution chain.
    print(f"'{old_english}' -> '{middle_english}' -> '{modern_english}'")
    print("-" * 50)

    final_word = modern_english
    print(f"Therefore, the Modern English word would likely be '{final_word}' (pronounced to rhyme with 'duster').")


if __name__ == '__main__':
    find_alternate_word()