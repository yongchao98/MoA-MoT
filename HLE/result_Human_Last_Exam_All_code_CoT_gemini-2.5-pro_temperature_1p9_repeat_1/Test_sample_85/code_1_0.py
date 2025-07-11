import sys
# This script explains the linguistic reconstruction of a Sanskrit word
# from its hypothetical Hindi and Punjabi descendants.

def explain_reconstruction():
    """
    Prints the step-by-step reasoning for the historical reconstruction.
    """
    hindi_word = "kãːs"
    punjabi_word = "kə̃d͡ʒʱ"
    
    print("Reconstructing the Sanskrit ancestor for Hindi '{}' and Punjabi '{}'.\n".format(hindi_word, punjabi_word))
    
    print("1. Initial Consonant and Vowel Analysis:")
    print("   - Both words start with 'k', so the ancestor also began with 'k'.")
    print("   - Both words have a nasalized vowel (Hindi ãː, Punjabi ə̃). This points to a Sanskrit structure of 'vowel + nasal consonant'.")
    print("   - The long vowel in Hindi (aː) versus the short vowel in Punjabi (ə) suggests compensatory lengthening in Hindi due to the simplification of a consonant cluster that Punjabi retained in a different form.\n")
    
    print("2. The Central Problem: Medial Consonant Divergence:")
    print("   - The most significant clue is the difference between Hindi '-s' and Punjabi '-d͡ʒʱ'.")
    print("   - We need a single Sanskrit consonant or cluster that could evolve into these two different sounds.\n")

    print("3. The Solution: The Sanskrit Cluster '-ṃś-'")
    print("   - The Sanskrit cluster '-ṃś-' (a nasal vowel followed by the palatal sibilant 'ś', IPA /ɕ/) provides the perfect source.")
    
    print("\n   - Path to Hindi:")
    print("     Skt. -ṃś- → Middle Indo-Aryan -ṃs- → Hindi -ãːs")
    print("     This is a regular sound change. The cluster simplifies to 's', and the preceding vowel becomes long and nasalized. Example: Skt. 'vaṃśa' → Hindi 'bā̃s'.")
    
    print("\n   - Path to Punjabi:")
    print("     Skt. -ṃś- → Punjabi -ə̃d͡ʒʱ-")
    print("     This change is more complex but also regular in Punjabi phonology. The sibilant 'ś' after a nasal is voiced and develops into an affricate with breathy voice ('d͡ʒʱ').\n")

    print("4. Final Reconstruction:")
    print("   - Combining the parts: k + a + ṃś + a (a common final vowel lost in the modern languages).")
    print("   - The reconstructed Sanskrit word is 'kaṃśa'.")
    
    # Python 2/3 compatibility for IPA characters
    if sys.version_info[0] < 3:
        reload(sys)
        sys.setdefaultencoding('utf-8')

    print("\n   - In the International Phonetic Alphabet (IPA), 'kaṃśa' is transcribed as: kãɕa")

explain_reconstruction()
