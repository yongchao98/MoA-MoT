def find_chomsky_nouns():
    """
    Identifies and prints the two nouns from the last syntactically correct and
    incorrect sentences in a key section of Chomsky's 'Syntactic Structures'.
    """
    last_correct_sentence = "sincerity may frighten the boy"
    last_incorrect_sentence = "sincerity may virtue the boy"

    noun1 = "sincerity"
    noun2 = "boy"

    print("In Noam Chomsky's 'Syntactic Structures' (1957), the chapter 'The Independence of Grammar' contains several example sentences.")
    print("-" * 80)

    print("The last syntactically CORRECT sentence he provides in the relevant list is:")
    print(f"  '{last_correct_sentence}'")
    print("\nThe last syntactically INCORRECT sentence is:")
    print(f"  '{last_incorrect_sentence}'")
    print("-" * 80)

    print("\nThe two nouns used in these sentences are:")
    print(f"1. {noun1}")
    print(f"2. {noun2}")

if __name__ == "__main__":
    find_chomsky_nouns()