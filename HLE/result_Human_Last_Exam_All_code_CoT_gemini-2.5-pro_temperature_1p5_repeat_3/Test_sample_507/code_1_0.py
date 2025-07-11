import pandas as pd

def find_similar_language_words():
    """
    This function presents a solution to the linguistic riddle by showcasing
    the surprising similarities between Korean and Tamil.
    """
    print("Finding two languages from different Asian cultures with no direct contact,")
    print("sharing similar words for 'mom', 'dad', and 'broom'.\n")
    print("The best candidates for this are Korean and Tamil (a Dravidian language from Southern India).")
    print("They belong to different language families and are geographically distant, yet share these similarities.\n")

    data = {
        'Concept': ['Language Family', 'Region', 'Mom', 'Dad', 'Broom (Root Word)'],
        'Korean': ['Koreanic (Isolate)', 'East Asia', 'eomma (엄마)', 'appa (아빠)', 'bi (비)'],
        'Tamil': ['Dravidian', 'South Asia', 'ammā (அம்மா)', 'appā (அப்பா)', 'pir (பிர்)*']
    }

    df = pd.DataFrame(data)

    print("--- Comparative Table ---")
    print(df.to_string(index=False))
    print("\n*Note: The similarity for 'broom' is based on a linguistic hypothesis connecting the Korean 'bi' (broom)")
    print("with the old Tamil root 'pir' (meaning to separate or divide), which is an act associated with sweeping.")
    print("The similarities for 'mom' and 'dad' are widely recognized and striking.")

    print("\n--- Word-by-Word Comparison ---")
    print(f"Mom: The Korean word '{data['Korean'][2]}' is nearly identical to the Tamil word '{data['Tamil'][2]}'.")
    print(f"Dad: The Korean word '{data['Korean'][3]}' is nearly identical to the Tamil word '{data['Tamil'][3]}'.")
    print(f"Broom: The Korean word '{data['Korean'][4]}' is considered similar to the Tamil root '{data['Tamil'][4]}' by some linguists.")

if __name__ == '__main__':
    find_similar_language_words()

<<<Korean and Tamil>>>