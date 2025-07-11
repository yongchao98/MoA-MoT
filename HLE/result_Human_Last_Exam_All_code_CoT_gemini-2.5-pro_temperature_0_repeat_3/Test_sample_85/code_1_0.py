import unicodedata

def reconstruct_sanskrit_ancestor():
    """
    This function reconstructs the Sanskrit ancestor for Hindi 'kā̃s' and Punjabi 'kə̃d͡ʒʱ'.
    
    The reconstruction is based on the following analysis:
    1.  Initial 'k' is consistent in both languages, pointing to Sanskrit 'k'.
    2.  The nasalized vowels (Hindi 'ãː', Punjabi 'ə̃') point to an original sequence of 'a' + Nasal consonant.
        The nasal would be the velar nasal 'ṅ' [ŋ] before a velar sound.
    3.  The core correspondence is Hindi '-s' vs. Punjabi '-d͡ʒʱ'. This points to the Sanskrit cluster '-kṣ-'.
        - Punjabi: The change Skt. -kṣ- > Pjb. -d͡ʒʱ- is attested (e.g., Skt. kakṣa > Pjb. kajjh).
        - Hindi: The change Skt. -kṣ- > H. -s is a plausible, though less common, development.
    4.  Combining these elements (k + aṅ + kṣ + a) gives the Sanskrit word 'kaṅkṣa'.
    """
    
    # The reconstructed Sanskrit word in standard academic transcription is "kaṅkṣa".
    # We will now represent this in the International Phonetic Alphabet (IPA).
    # k = Voiceless velar plosive
    # a = Open central unrounded vowel (often transcribed as 'ɐ')
    # ŋ = Voiced velar nasal
    # k = Voiceless velar plosive
    # ʂ = Voiceless retroflex fricative
    # a = Open central unrounded vowel (schwa, often transcribed as 'ɐ')
    
    reconstructed_ipa = "kaŋkʂa"
    
    # The problem asks to output each number in the final equation, which seems to be a misplaced instruction.
    # I will interpret it as printing the components of the reconstruction.
    print("Reconstruction Steps:")
    print("Initial Consonant: k")
    print("Vowel + Nasal: aŋ (from a + ṅ)")
    print("Medial Cluster: kʂ")
    print("Final Vowel: a")
    print("\nFinal Reconstructed Word (IPA):")
    print(reconstructed_ipa)

reconstruct_sanskrit_ancestor()