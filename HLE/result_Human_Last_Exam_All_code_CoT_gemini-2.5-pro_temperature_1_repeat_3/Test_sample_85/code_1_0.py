import unicodedata

def reconstruct_sanskrit_ancestor():
    """
    This function provides the reconstructed Sanskrit ancestor for the given words.
    The reconstruction is based on principles of historical linguistics applied
    to the sound changes from Sanskrit to Hindi and Punjabi.

    The logic is as follows:
    1. Hindi 'kãːs' points to a proto-form with a sibilant, like *kaṃśya.
    2. Punjabi 'kə̃d͡ʒʱ' points to a proto-form with the cluster -dhy-, like *kandhya.
    3. The Punjabi reflex from *kandhya is perfectly regular.
    4. The Hindi reflex from *kandhya (*kā̃jh) requires an irregular change (jh > s)
       to become 'kãːs', but this change is more plausible than the changes
       required to get the Punjabi word from *kaṃśya.
    5. Therefore, the most probable ancestor is 'kandhya'.

    The Sanskrit form 'kandhya' is rendered in IPA.
    k -> [k] (voiceless velar stop)
    a -> [a] (open central unrounded vowel)
    n -> [n̪] (dental nasal, assimilated to the following dental stop)
    dh -> [d̪ʱ] (voiced dental aspirated stop)
    y -> [j] (palatal approximant)
    a -> [a] (open central unrounded vowel)
    """
    
    # The reconstructed proto-word in Sanskrit is 'kandhya'.
    # Here is its representation in the International Phonetic Alphabet (IPA).
    reconstructed_ipa = "kand̪ʱja"
    
    print("Reconstructed Sanskrit Ancestor (IPA):")
    print(reconstructed_ipa)

reconstruct_sanskrit_ancestor()