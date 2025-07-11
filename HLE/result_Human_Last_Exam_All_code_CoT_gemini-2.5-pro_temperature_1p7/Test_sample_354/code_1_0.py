import textwrap

def solve_disneyfication_question():
    """
    Explains the dimensions of Disneyization and determines the correct answer.
    """
    # Alan Bryman's theory of Disneyization (or Disneyfication) outlines four
    # key processes by which the principles of Disney theme parks are applied
    # to other sectors of society.

    dimensions = {
        "Theming": "The application of a narrative to an institution or location to create a holistic and immersive experience.",
        "Hybrid Consumption": "The blurring of lines between different forms of consumption, where, for example, a store is also an entertainment venue and a restaurant.",
        "Merchandising": "The promotion and sale of goods bearing promotional images or logos related to the theme or brand.",
        "Performative Labor": "The scripting of employees' words and actions, and the management of their emotional expressions to enhance the customer experience. Employees become 'cast members' in a performance."
    }

    print("Alan Bryman's four dimensions of Disneyization are:\n")
    for dim, desc in dimensions.items():
        print(f"- {dim}:")
        print(textwrap.fill(desc, width=80, initial_indent='  ', subsequent_indent='  '))
        print()

    # Let's analyze the given options:
    # A. hybrid consumption and merchandising - Both are correct dimensions.
    # B. performative labor and sanitization - 'Sanitization' is a result of the process, not a core dimension.
    # C. trivialization and theming - 'Trivialization' is a common critique, not a core dimension.
    # D. sanitization and trivialization - Neither are core dimensions.
    # E. Disneyfication and Disneyization - These are just names for the concept.
    # F. McDonaldization and Disneyization - These are two different sociological concepts.
    # G. theming and performative labor - Both are correct dimensions.

    # Both options A and G list two correct dimensions from Bryman's theory.
    # However, 'theming' and 'performative labor' are often considered the most
    # foundational aspects that deal with the transformation of social space
    # and human interaction, which is at the heart of the "Disneyization of Society."

    correct_answer = "G"

    print("--------------------------------------------------")
    print(f"The most accurate choice that represents two core characteristics of Disneyfication discussed by Bryman is G: theming and performative labor.")

solve_disneyfication_question()