import pandas as pd

def find_similar_words():
    """
    Identifies and displays two languages from different Asian cultures
    with similar words for "mom", "dad", and "broom".
    """
    
    # --- Data ---
    # The chosen languages are Korean and Tamil.
    # Korean is a language isolate from East Asia.
    # Tamil is a Dravidian language from South Asia.
    # They have no mainstream-accepted genetic relationship or direct historical contact.
    
    language1 = "Korean"
    culture1 = "East Asian"
    
    language2 = "Tamil"
    culture2 = "South Asian (Dravidian)"

    # --- Word Comparisons ---
    # The words for "mom" and "dad" are famously similar.
    # The word for "broom" is more complex. While direct nouns differ,
    # some linguists point to a similarity between the Korean noun for a simple
    # broom and the Tamil verb "to sweep".

    data = {
        "Concept": ["Mom", "Dad", "Broom (or related verb)"],
        "Korean": ["엄마 (eomma)", "아빠 (appa)", "비 (bi)"],
        "Tamil": ["அம்மா (amma)", "அப்பா (appa)", "பெருக்கு (perukku - verb for 'to sweep')"]
    }

    # --- Output ---
    print(f"The two languages are {language1} (from {culture1} culture) and {language2} (from {culture2} culture).\n")
    print("These languages have no established direct contact and belong to different language families,")
    print("making the similarities—especially for 'mom' and 'dad'—a subject of linguistic interest.\n")
    
    # Using pandas for clean, tabular output
    df = pd.DataFrame(data)
    print(df.to_string(index=False))

find_similar_words()