import unicodedata

def reconstruct_sanskrit_ancestor():
    """
    This function reconstructs the Sanskrit ancestor for Hindi kãːs and Punjabi kə̃d͡ʒʱ.

    The reconstruction hinges on the divergent evolution of a specific Sanskrit consonant cluster.

    1.  k- -> k-: The initial consonant is stable in both languages.
    2.  -ãː- / -ə̃-: The nasalized vowel in both descendants points to an original short vowel 'a' followed by a cluster containing a nasal, like *aNC. The vowel length difference is a regular result of different cluster simplification rules (compensatory lengthening in Hindi, shortening in Punjabi).
    3.  -s- / -d͡ʒʱ-: This is the key correspondence. This points to the versatile Sanskrit cluster *-ṅkṣ- (/ŋkʂ/).
        - Path to Hindi: Skt. *-ṅkṣ- > Prakrit *-ṃcch- > Hindi *-̃ch- > Hindi *-̃s-. The change from 'ch' to 's' is a known simplification. e.g., pūchnā > pūsnā ('to ask').
        - Path to Punjabi: Skt. *-ṅkṣ- > *-gž- > *-dʒ-. The voicing of 'kṣ' is a known feature in some Northwestern Indo-Aryan languages. The aspiration in 'd͡ʒʱ' is likely a secondary, spontaneous development.

    Combining these elements, the most probable ancestor is 'kaṅkṣa'.
    """

    # The reconstructed word in Devanagari is कङ्क्ष
    # The standard academic transcription is 'kaṅkṣa'
    reconstructed_ipa = "kaŋkʂa"

    print("Reconstructed Sanskrit Ancestor (in IPA):")
    print(reconstructed_ipa)

reconstruct_sanskrit_ancestor()