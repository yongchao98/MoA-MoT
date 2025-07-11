import textwrap

def solve_linguistic_puzzle():
    """
    Analyzes a list of words to identify which one has not undergone
    Trisyllabic Laxing (TSL) in its sound development.
    The analysis and conclusion are printed to the console.
    """

    # Step 1: Define the core concept
    print("--- Analyzing Words for Trisyllabic Laxing (TSL) ---")
    print("\n")
    
    tsl_definition = (
        "Trisyllabic Laxing is a rule in English phonology where a tense vowel "
        "(like /aɪ/ in 'divine') becomes a lax vowel (like /ɪ/ in 'divinity') "
        "when it is in the antepenultimate (third-from-last) syllable of a word."
    )
    print("Definition:")
    print(textwrap.fill(tsl_definition, width=80))
    print("\n" + "="*80 + "\n")

    # Step 2: Create a data structure for analysis
    # We will analyze each word against its base form to check for vowel laxing.
    word_analysis = {
        'derivative': "HAS undergone TSL. The base word is 'derive' with a tense vowel /aɪ/. In 'derivative' /dɪ-'rɪv-ə-tɪv/, the vowel laxes to /ɪ/ because it is followed by two syllables.",
        'serenity': "HAS undergone TSL. The base word is 'serene' with a tense vowel /iː/. In 'serenity' /sə-'rɛn-ɪ-ti/, the vowel laxes to /ɛ/ because it is followed by two syllables.",
        'pleasant': "HAS undergone a TSL-like process. The base word is 'please' with a tense vowel /iː/. In 'pleasant' /'plɛz-ənt/, the vowel laxes to /ɛ/. This is a common pattern of shortening.",
        'gratitude': "HAS undergone TSL. We compare it to the related word 'grateful,' which shows the tense vowel /eɪ/. In 'gratitude' /'græt-ɪ-tud/, the vowel laxes to /æ/ before the following two syllables.",
        'southern': "has NOT undergone Trisyllabic Laxing. While the vowel changes from the base 'south' (/aʊ/ -> /ʌ/), it is only followed by ONE syllable (-ern). The rule requires two following syllables, so this is a different sound change, not TSL.",
        'shadow': "has NOT undergone TSL. This word comes from Old English 'sceadu,' which already had a short vowel /æ/. Since there was no tense vowel in its historical source to begin with, the 'laxing' process could not have occurred."
    }

    print("Step-by-Step Analysis:\n")

    # Step 3: Print the analysis and conclusion
    # The final equation here is the logical deduction process.
    answer = None
    for word, analysis in word_analysis.items():
        print(f"Word: {word}")
        print(f"  -> Analysis: {analysis}\n")
        # 'shadow' is the most definitive example of not undergoing the process
        if "NOT undergone TSL" in analysis and word == "shadow":
            answer = word

    print("="*80 + "\n")
    print("Conclusion:")
    print("While 'southern' doesn't fit the syllable count for TSL, 'shadow' is a clearer example of a word that was never eligible for the process at all.")
    print("Its vowel was already short in its Old English source word, so it never had a tense vowel that could be laxed.")
    
    print(f"\nFinal Answer Equation: shadow (source vowel was already lax) ≠ derivative (tense -> lax)")
    
solve_linguistic_puzzle()
<<<shadow>>>