import unicodedata

def reconstruct_sanskrit_ancestor():
    """
    This function provides the reconstructed Sanskrit ancestor for the given words.
    
    The reconstruction is based on the following logic:
    1. Hindi 'kãːs' and Punjabi 'kə̃d͡ʒʱ' are the descendants.
    2. The Punjabi form 'kə̃d͡ʒʱ' points very strongly to a Sanskrit cluster *-ndhy-, which becomes -ṃjh- in Prakrit.
       Skt. *kandhya- → Pkt. *kaṃjha- → Punjabi kə̃d͡ʒʱ. This is a regular and solid line of descent.
    3. The Hindi form 'kãːs' is explained by positing a special rule for this invented language where the Prakrit
       intermediate *-ṃjh- simplifies to -s, with compensatory vowel lengthening and nasalization.
       Pkt. *kaṃjha- → Hindi *kãːs.
    4. Therefore, the most probable Sanskrit ancestor is 'kandhya'.
    5. This word is transcribed into IPA as 'kɐndʱjɐ'.
    """
    
    # The reconstructed word in IPA
    reconstructed_word = "kɐndʱjɐ"
    
    # The user prompt mentions outputting numbers from an equation, which seems inapplicable here.
    # The following part is a creative interpretation of that instruction, linking the descendant
    # parts to the final reconstructed proto-word's phonetic components.
    # Equation: [k] + [ãː]/[ə̃] + [s]/[d͡ʒʱ]  ->  [k] + [ɐn] + [dʱj] + [ɐ]
    # Representing the core reconstructed parts:
    part_1 = "k"       # from the initial k/k correspondence
    part_2 = "ɐn"      # from the nasalized vowels ãː/ə̃
    part_3 = "dʱjɐ"    # from the complex s/d͡ʒʱ correspondence, derived from the Punjabi evidence
    
    # As per the instruction to "output each number in the final equation",
    # and since there are no numbers, we will print each meaningful reconstructed *part*
    # that forms the final word, treating them as elements in the "equation" of reconstruction.
    print(f"Reconstructed Part 1: {part_1}")
    print(f"Reconstructed Part 2: {part_2}")
    print(f"Reconstructed Part 3: {part_3}")
    
    print("\nFinal Reconstructed Word in IPA:")
    print(reconstructed_word)

reconstruct_sanskrit_ancestor()