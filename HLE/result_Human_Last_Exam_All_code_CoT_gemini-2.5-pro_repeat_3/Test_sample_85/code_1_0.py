# This script reconstructs the Sanskrit ancestor of two given Indo-Aryan words.

def reconstruct_sanskrit_ancestor():
    """
    Performs the linguistic reconstruction based on phonological rules
    and prints the final result in IPA.
    """

    # Input words from two descendant languages.
    hindi_word = "kãːs"
    punjabi_word = "kə̃d͡ʒʱ"

    # Step 1: Analyze the initial consonant.
    # 'k' is stable in both Hindi and Punjabi, corresponding to Sanskrit 'k'.
    reconstructed_initial = "k"

    # Step 2: Analyze the vowel and nasalization.
    # Hindi 'ãː' and Punjabi 'ə̃' point to an original short vowel 'a'
    # followed by a nasal consonant.
    reconstructed_vowel = "a"
    # The nasal is part of the medial cluster.

    # Step 3: Analyze the medial consonant correspondence.
    # The correspondence is Hindi '-s-' vs. Punjabi '-d͡ʒʱ-'.
    # The Punjabi voiced aspirate '-d͡ʒʱ-' combined with nasalization strongly
    # suggests an original cluster of 'Nasal + Voiceless Consonant'.
    # The Sanskrit cluster '-kṣ-' becomes '-cch-' in Punjabi.
    # A preceding nasal voices this cluster: N + cch -> ñjh [d͡ʒʱ].
    # Therefore, the source cluster is most likely '-ṅkṣ-'.
    # In some Hindi dialects, '-kṣ-' simplifies to '-s-'.
    # Thus, 'kaṅkṣa' > '*kaṁsa' > 'kãːs' is a valid path for Hindi.
    reconstructed_cluster = "ṅkṣ" # Using Devanagari representation for clarity

    # Step 4: Assemble the final reconstruction.
    # The full form is ka + ṅkṣ + a.
    sanskrit_reconstruction = "kaṅkṣa"

    # Step 5: Convert the final reconstruction to International Phonetic Alphabet (IPA).
    # k -> k
    # a -> a
    # ṅ -> ŋ (velar nasal)
    # kṣ -> kʂ (cluster of 'k' and retroflex sibilant 'ʂ')
    # a -> a
    ipa_reconstruction = "kaŋkʂa"

    # Print the final reconstructed ancestor in IPA.
    # The prompt asks to "output each number in the final equation", which is
    # likely a misplaced instruction. We will print the final string clearly.
    print("Reconstructed Sanskrit Ancestor (in IPA):")
    print(ipa_reconstruction)

reconstruct_sanskrit_ancestor()
<<<kaŋkʂa>>>