import sys

# In Python 3, to ensure unicode characters like 'ā' or 'ī' print correctly on all systems.
if sys.stdout.encoding != 'UTF-8':
    sys.stdout.reconfigure(encoding='utf-8')

def analyze_words():
    """
    Analyzes a list of words to determine which has not undergone trisyllabic laxing.
    """
    words_data = {
        'southern': {'syllables': 2, 'root': 'south', 'root_vowel': 'aʊ (as in mouth)', 'word_vowel': 'ʌ (as in strut)', 'origin': 'Germanic'},
        'derivative': {'syllables': 4, 'root': 'derive', 'root_vowel': 'aɪ (as in prize)', 'word_vowel': 'ɪ (as in kit)', 'origin': 'Latinate'},
        'serenity': {'syllables': 4, 'root': 'serene', 'root_vowel': 'iː (as in fleet)', 'word_vowel': 'ɛ (as in dress)', 'origin': 'Latinate'},
        'pleasant': {'syllables': 2, 'root': 'please', 'root_vowel': 'iː (as in fleet)', 'word_vowel': 'ɛ (as in dress)', 'origin': 'Latinate'},
        'gratitude': {'syllables': 3, 'root': 'grateful/grace', 'root_vowel': 'eɪ (as in face)', 'word_vowel': 'æ (as in trap)', 'origin': 'Latinate'},
        'shadow': {'syllables': 2, 'root': 'shade', 'root_vowel': 'eɪ (as in face)', 'word_vowel': 'æ (as in trap)', 'origin': 'Germanic'}
    }

    print("Analyzing each word for Trisyllabic Laxing (TSL)...")
    print("Rule: TSL shortens a long vowel in the 3rd-to-last syllable of a word (word must have >= 3 syllables).\n")
    
    non_tsl_candidates = []

    for word, data in words_data.items():
        print(f"--- Analyzing: {word.capitalize()} ---")
        print(f"Syllable count: {data['syllables']}")
        print(f"Comparison: Root '{data['root']}' [{data['root_vowel']}] -> '{word}' [{data['word_vowel']}]")
        
        if data['syllables'] >= 3:
            # All words with >= 3 syllables in this list are textbook examples of TSL.
            # The laxing happens in the antepenultimate syllable.
            # de-RIV-a-tive, se-REN-i-ty, GRAT-i-tude
            print("Verdict: MEETS the syllable requirement for TSL. Vowel shortening occurs in the antepenultimate syllable.")
            print("This word HAS undergone trisyllabic laxing.")
        else:
            print(f"Verdict: DOES NOT MEET the syllable requirement for TSL (has only {data['syllables']} syllables).")
            print("This word has NOT undergone trisyllabic laxing.")
            non_tsl_candidates.append(word)
        print("-" * (len(word) + 20) + "\n")

    print("--- Conclusion ---")
    print(f"The words that have not undergone trisyllabic laxing are: {', '.join(non_tsl_candidates)}.")
    print("All three fail because they have fewer than three syllables.")
    print("However, since a single answer is required, a distinction can be made.")
    print("Trisyllabic laxing is a process that primarily applies to words of Latin/Romance origin.")
    print(f"'{words_data['derivative']['origin']}', '{words_data['serenity']['origin']}', '{words_data['gratitude']['origin']}', and '{words_data['pleasant']['origin']}' are all Latinate.")
    print(f"'southern' and 'shadow' are Germanic, meaning their vowel changes stem from different historical processes.")
    print("Of these, 'southern' is a classic linguistic example of vowel shortening that contrasts with 'northern' (which does not shorten), making it a strong example of a non-TSL process.")
    final_answer = "southern"
    print(f"\nThe best answer is '{final_answer}'.")
    
if __name__ == '__main__':
    analyze_words()
    print("\n<<<southern>>>")