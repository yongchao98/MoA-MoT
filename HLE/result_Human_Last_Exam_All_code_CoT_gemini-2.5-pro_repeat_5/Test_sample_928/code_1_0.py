def find_chomsky_nouns():
    """
    This script identifies the two nouns from the last syntactically correct and
    incorrect sentences in the relevant section of Chomsky's "Syntactic Structures".
    """

    # The last pair of sentences used by Chomsky in the section to
    # demonstrate a grammatical/ungrammatical distinction.
    last_correct_sentence = "I saw a fragile whale"
    last_incorrect_sentence = "I saw a fragile of"

    # Identify nouns (including pronouns) in the correct sentence.
    # "I" is a pronoun, "whale" is a common noun.
    nouns_in_correct = ["I", "whale"]

    # Identify nouns (including pronouns) in the incorrect sentence.
    # "I" is a pronoun. The sentence is ungrammatical due to the lack of a noun after "of".
    nouns_in_incorrect = ["I"]

    # The question asks for the two nouns used in these sentences.
    # Combining the unique nouns from both gives us the answer.
    unique_nouns = set(nouns_in_correct + nouns_in_incorrect)
    
    noun1, noun2 = unique_nouns

    print(f"The last syntactically correct sentence Chomsky gives is: '{last_correct_sentence}'")
    print(f"The noun and pronoun used in this sentence are: '{nouns_in_correct[0]}' and '{nouns_in_correct[1]}'")
    print("")
    print(f"The last syntactically incorrect sentence Chomsky gives is: '{last_incorrect_sentence}'")
    print(f"The pronoun used in this sentence is: '{nouns_in_incorrect[0]}'")
    print("")
    print(f"The two distinct nouns (including pronouns) used across both sentences are:")
    print(f"1. {noun1}")
    print(f"2. {noun2}")

find_chomsky_nouns()