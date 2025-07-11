import pandas as pd

def analyze_trisyllabic_laxing():
    """
    Analyzes a list of words to determine which has not undergone Trisyllabic Laxing.
    Trisyllabic Laxing is a sound change where a tense vowel becomes lax when it is
    in the third-to-last syllable of a word.
    """
    
    data = {
        'Word': ['southern', 'derivative', 'serenity', 'pleasant', 'gratitude', 'shadow'],
        'Base/Root': ['south', 'derive', 'serene', 'please', 'grate', "Old English 'sceadu'"],
        'Base Vowel': ['/aʊ/ (tense)', '/aɪ/ (tense)', '/iː/ (tense)', '/iː/ (tense)', '/eɪ/ (tense)', '/æ/ (lax)'],
        'Derived Vowel': ['/ʌ/ (lax)', '/ɪ/ (lax)', '/ɛ/ (lax)', '/ɛ/ (lax)', '/æ/ (lax)', '/æ/ (lax)'],
        'Underwent Laxing?': ['Yes (shortening)', 'Yes (TSL)', 'Yes (TSL)', 'Yes (shortening)', 'Yes (TSL)', 'No']
    }

    analysis_df = pd.DataFrame(data)
    
    print("--- Analysis of Vowel Laxing in Each Word ---")
    
    final_answer = ""
    for index, row in analysis_df.iterrows():
        word = row['Word']
        base_vowel = row['Base Vowel']
        derived_vowel = row['Derived Vowel']
        underwent_laxing = row['Underwent Laxing?']
        
        # We want to display the final equation for each word.
        # This will be in the format: "Base Vowel -> Derived Vowel"
        base_vowel_sound = base_vowel.split(' ')[0]
        derived_vowel_sound = derived_vowel.split(' ')[0]
        
        print(f"\nWord: {word.upper()}")
        
        if underwent_laxing == 'No':
            print(f"  Analysis: The root vowel {base_vowel_sound} was already lax and did not change.")
            print(f"  Final Equation: {base_vowel_sound} -> {derived_vowel_sound} (no change from tense to lax)")
            final_answer = word
        else:
            print(f"  Analysis: The original tense vowel {base_vowel_sound} became the lax vowel {derived_vowel_sound}.")
            print(f"  Final Equation: {base_vowel_sound} -> {derived_vowel_sound}")

    print("\n--- Conclusion ---")
    print(f"All words except '{final_answer}' show a historical tense-to-lax vowel shift.")
    print(f"The word '{final_answer}' did not undergo this process because its vowel was historically already lax.")

    return final_answer

if __name__ == '__main__':
    result = analyze_trisyllabic_laxing()
    print(f"\n<<<{result}>>>")