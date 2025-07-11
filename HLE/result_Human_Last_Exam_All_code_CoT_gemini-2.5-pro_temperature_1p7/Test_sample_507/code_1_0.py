import pandas as pd

def find_linguistic_coincidence():
    """
    This function explains and presents the solution to the linguistic puzzle.
    It identifies two languages from two different Asian cultures with no direct
    contact that share similar words for "mom", "dad", and "broom".
    """

    print("Searching for two languages in Asia from different cultures with no direct contact...")
    print("The words for 'mom', 'dad', and 'broom' must be highly similar.\n")

    # The chosen languages
    lang1 = "Korean"
    culture1 = "East Asia"
    lang2 = "Persian (Farsi)"
    culture2 = "West Asia"

    # Explanation of the choice
    print(f"Candidate Pair Found: {lang1} ({culture1}) and {lang2} ({culture2})")
    print("These languages belong to different language families (Koreanic and Indo-European)")
    print("and their cultures developed thousands of miles apart with no direct contact.\n")

    # Create a DataFrame for clear comparison
    data = {
        "Word": ["Mom", "Dad", "Broom"],
        "Korean": ["eomma (엄마)", "appa (아빠)", "bitjaru (빗자루)"],
        "Persian (Farsi)": ["maman (مامان)", "baba (بابا)", "jâru (جارو)"]
    }
    df = pd.DataFrame(data)

    print("Comparing the words:\n")
    print(df.to_string(index=False))
    print("\n")

    # Detailed explanation of similarities
    print("Analysis of Similarity:")
    print("-----------------------")
    print("'Mom'  : 'eomma' and 'maman' share the 'm' sound, common in words for mother globally.")
    print("'Dad'  : 'appa' and 'baba' both use bilabial stops ('p' or 'b'), a common feature in 'papa/dada'-type words.")
    print("'Broom': The similarity is striking here. The Persian word is 'jâru', and the Korean word 'bitjaru' contains an identical-sounding component, 'jaru'. The Korean word breaks down into 'bit' (comb) + 'jaru' (handle). The root words in each language are unrelated, making this a remarkable coincidence.")

    # Final answer format as requested
    final_answer = f"{lang1} and {lang2}"
    print(f"\n<<<The two languages are {final_answer}>>>")

if __name__ == '__main__':
    find_linguistic_coincidence()