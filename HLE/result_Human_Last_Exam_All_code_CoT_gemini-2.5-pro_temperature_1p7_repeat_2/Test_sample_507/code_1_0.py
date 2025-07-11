import sys

def find_similar_words():
    """
    This function presents the findings for two Asian languages with similar words
    for "broom", "mom", and "dad" despite having no direct contact.
    """
    
    lang1 = "Korean"
    lang2 = "Persian (Farsi)"

    print(f"The two languages are {lang1} and {lang2}.")
    print("\nThese languages are from distinct cultures in East and West Asia and belong to different language families.")
    print("The similarities are considered coincidental or due to universal linguistic patterns (especially for 'mom' and 'dad').")
    print("-" * 40)
    
    # --- Word Data ---
    words = {
        "Broom": {
            "lang1_word": "(bit)jaru",
            "lang1_desc": "빗자루 (bit-jaru)",
            "lang2_word": "jâru",
            "lang2_desc": "جارو (jâru)",
            "note": "The Korean word is a compound; the 'jaru' part, meaning 'handle', is nearly identical to the Persian word."
        },
        "Mom": {
            "lang1_word": "eomma",
            "lang1_desc": "엄마 (eomma)",
            "lang2_word": "māmān",
            "lang2_desc": "مامان (māmān)",
            "note": "Both are 'nursery words' based on the bilabial nasal sound 'm'."
        },
        "Dad": {
            "lang1_word": "appa",
            "lang1_desc": "아빠 (appa)",
            "lang2_word": "bābā",
            "lang2_desc": "بابا (bābā)",
            "note": "Both are 'nursery words' based on the repetitive bilabial plosive sound 'p' or 'b'."
        }
    }

    # --- Print the Comparisons ---
    for i, (english_word, data) in enumerate(words.items(), 1):
        print(f"\n{i}. Word: {english_word}")
        print(f"   {lang1}: {data['lang1_desc']}")
        print(f"   {lang2}: {data['lang2_desc']}")
        print("\n   Comparison of the core sounds:")
        print(f"   {data['lang1_word']}  <-->  {data['lang2_word']}")
        print(f"   Note: {data['note']}")

find_similar_words()

# Appending the final answer in the specified format to stdout.
# This ensures it's the last part of the output.
sys.stdout.write("\n<<<Korean and Persian (Farsi)>>>")