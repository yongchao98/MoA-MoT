import textwrap

# --- Introduction ---
language1 = "Turkish"
language2 = "Malay"

intro_text = f"""
Two languages from different Asian cultures with no direct linguistic link, yet sharing remarkably similar words for "mom," "dad," and "broom," are {language1} and {language2}.

{language1} is part of the Turkic language family, primarily spoken in Western Asia. {language2} belongs to the Austronesian family, spoken in Maritime Southeast Asia. Their similarities are a mix of universal linguistic patterns and pure coincidence.
"""
print(intro_text)
print("-" * 50)

# --- Word Data & Comparison ---
comparisons = [
    {
        "word": "DAD",
        "lang1_word": "baba",
        "lang2_word": "bapak",
        "analysis": "These words are nearly identical. This is a classic example of a 'nursery word' or 'Lallwort', based on the simple labial consonant 'b' and vowel 'a', which are among the first sounds babies can make."
    },
    {
        "word": "MOM",
        "lang1_word": "ana",
        "lang2_word": "emak",
        "analysis": "While not identical, both 'ana' (the Turkic root for mother) and 'emak' (a common word for mother in Malay) are built on a similar 'nursery word' principle, using a nasal consonant ('n' or 'm') and an open vowel."
    },
    {
        "word": "BROOM",
        "lang1_word": "süpürge",
        "lang2_word": "sapu",
        "analysis": "This is the most striking similarity because it's purely coincidental. The Turkish root 'süpür-' (to sweep) and the Malay word 'sapu' (to sweep/broom) sound very similar but developed independently from their respective proto-languages (Proto-Turkic and Proto-Malayo-Polynesian)."
    }
]

for item in comparisons:
    print(f"\nComparing the word for '{item['word']}':")
    print(f"  - In {language1}:  '{item['lang1_word']}'")
    print(f"  - In {language2}:  '{item['lang2_word']}'")
    print("\n  Analysis:")
    # Wrap the analysis text for better readability
    wrapped_analysis = textwrap.fill(item['analysis'], width=60, initial_indent='    ', subsequent_indent='    ')
    print(wrapped_analysis)
    print("-" * 50)
