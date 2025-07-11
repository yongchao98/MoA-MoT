def solve_reconstruction():
    """
    Reconstructs the Sanskrit ancestor for Hindi kãːs and Punjabi kə̃dʒʱ
    and explains the step-by-step reasoning.
    """
    
    hindi_word = "kãːs"
    punjabi_word = "kə̃dʒʱ"
    reconstructed_ipa = "kanɕa"
    
    print("Reconstruction of the Sanskrit Ancestor")
    print("-" * 40)
    print(f"Descendant 1 (Hindi): {hindi_word}")
    print(f"Descendant 2 (Punjabi): {punjabi_word}\n")
    
    print("Step 1: Analyze the initial consonant.")
    print("The initial sound 'k' is identical in both Hindi (k) and Punjabi (k).")
    print("Therefore, the original Sanskrit word also began with 'k'.")
    print("Proto-form so far: k...")
    print("-" * 40)
    
    print("Step 2: Analyze the vowel and nasalization.")
    print("Hindi has a long nasal vowel ('ãː') while Punjabi has a short nasal vowel ('ə̃').")
    print("This pattern typically points to a Sanskrit short vowel ('a') followed by a nasal consonant cluster.")
    print("In the path to Hindi, the cluster was simplified, and the vowel was lengthened to compensate ('compensatory lengthening').")
    print("In the path to Punjabi, the vowel remained short.")
    print("Proto-form so far: ka-cluster-a")
    print("-" * 40)

    print("Step 3: Analyze the medial consonant correspondence.")
    print("This is the central clue. Hindi has 's' where Punjabi has 'd͡ʒʱ'.")
    print("We need to find a single Sanskrit cluster that could split and result in these two different sounds.")
    print("-" * 40)
    
    print("Step 4: Propose the proto-cluster '-nś-'.")
    print("The cluster '-nś-' (a nasal 'n' plus a palatal sibilant 'ś') provides the best solution.")
    print("-" * 40)
    
    print("Step 5: Trace the development to Hindi.")
    print("Sanskrit '-anśa' > Middle Indo-Aryan '-aṁsa' > Modern Hindi '-ãːs'")
    print("The '-nś-' cluster regularly simplifies to '-s-' in Hindi, with the nasalization moving to the vowel and lengthening it.")
    print(f"This perfectly explains the Hindi form: k + ãː + s = {hindi_word}.")
    print("-" * 40)

    print("Step 6: Trace the development to Punjabi.")
    print("Sanskrit '-anśa' > Middle Indo-Aryan '-aṁjha' (a rarer sound change) > Modern Punjabi '-ə̃d͡ʒʱ'")
    print("In the dialect leading to Punjabi, the palatal cluster '-nś-' was voiced and became an aspirated affricate '-jh-' (IPA d͡ʒʱ). The preceding vowel remained short.")
    print(f"This perfectly explains the Punjabi form: k + ə̃ + d͡ʒʱ = {punjabi_word}.")
    print("-" * 40)

    print("Step 7: Final Reconstruction.")
    print("Combining all elements, the reconstructed Sanskrit word is 'kanśa'.")
    print("In the International Phonetic Alphabet (IPA), 'ś' is written as 'ɕ'.")
    print("The final reconstructed form is represented in IPA as follows:")
    
    final_equation = list(reconstructed_ipa)
    print(" + ".join(final_equation))
    
    print(f"\nFinal Answer (as a single IPA string):")


solve_reconstruction()

# The final answer is wrapped in <<<>>>
print(f'<<<kanɕa>>>')
