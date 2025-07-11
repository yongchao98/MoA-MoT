def find_chomsky_nouns():
    """
    This function identifies and prints the two nouns from the last syntactically
    correct and last syntactically incorrect sentences in a key section of
    Noam Chomsky's "Syntactic Structures".

    The sentences are:
    - Last correct: "the boy may admire sincerity"
    - Last incorrect: "*the boy may virtue sincerity"
    """

    # Nouns identified from the sentences
    noun1 = "sincerity"
    noun2 = "boy"

    print("The two nouns from the specified Chomsky sentences are:")
    print(f"1. {noun1}")
    print(f"2. {noun2}")

if __name__ == "__main__":
    find_chomsky_nouns()