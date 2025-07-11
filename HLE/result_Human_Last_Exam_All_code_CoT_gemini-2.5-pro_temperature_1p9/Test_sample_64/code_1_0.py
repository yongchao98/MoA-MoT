import sys

def find_pitch_accent(word, reading, morae_count, accent_kernel):
    """
    Determines the Japanese pitch accent pattern type based on word properties.
    """
    # Define accent pattern rules based on morae count and kernel position
    if accent_kernel == 0:
        pattern_name = "Heiban"
    elif accent_kernel == 1:
        pattern_name = "Atamadaka"
    elif accent_kernel == morae_count:
        pattern_name = "Odaka"
    elif 1 < accent_kernel < morae_count:
        pattern_name = "Nakadaka"
    else:
        pattern_name = "Unknown"
        
    # --- Explanation ---
    print(f"Analyzing the pitch accent for the word: 「{word}」 ({reading})")
    print("-" * 30)
    print(f"1. The word has {morae_count} morae (o-to-u-to).")
    print(f"2. The standard accent kernel position is [{accent_kernel}].")
    print("3. This means the pitch rises after the first mora and falls after the third.")
    print("   Pitch visualization: Low-High-High-Low (おとうと)")
    print("-" * 30)
    print("Classifying the pattern:")
    print(f" - For a {morae_count}-mora word, an accent at [1] would be Atamadaka.")
    print(f" - An accent at [{morae_count-1}] or earlier (but not [1]) is Nakadaka.")
    print(f" - An accent at [{morae_count}] would be Odaka.")
    print("-" * 30)
    
    # --- Final Conclusion ---
    print("Final Equation:")
    print(f"Word morae = {morae_count}, Accent kernel = {accent_kernel} => The pattern is {pattern_name}.")
    
    return pattern_name

# Word information for 「弟」
word_kanji = "弟"
word_reading = "おとうと"
word_morae_count = 4
word_accent_kernel = 3 # Accent on the 3rd mora (とう)

# Run the analysis
find_pitch_accent(word_kanji, word_reading, word_morae_count, word_accent_kernel)
