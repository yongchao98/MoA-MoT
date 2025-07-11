def solve_solaris_question():
    """
    This script determines which character in the 1972 movie 'Solaris'
    is ashamed to miss the sound of leaves rustling on Earth.
    """
    
    question = "In the 1972 Andrei Tarkovsky movie Solaris which character is ashamed to miss the sound of leaves rustling on Earth?"
    
    options = {
        'A': 'Kris',
        'B': 'Hari',
        'C': 'Snaut',
        'D': 'Sartorius',
        'E': 'Gibarian'
    }

    # Step 1: Provide the detailed reasoning for the correct answer.
    explanation = (
        "In the film, the character Dr. Snaut has a conversation with Kris Kelvin where he discusses\n"
        "the psychological toll their long mission has taken. He has become weary and philosophical.\n"
        "In a moment of vulnerability, Snaut confesses his homesickness and specifically states\n"
        "that he is ashamed to admit he misses the sound of rustling leaves on Earth.\n"
        "This highlights a central theme: the crew's disconnection from nature and their own humanity."
    )

    # Step 2: Identify the correct option based on the reasoning.
    correct_option = 'C'
    correct_character = options[correct_option]

    # Step 3: Print the analysis and the final conclusion.
    # The format below emulates an "equation" as requested, showing the logic flow.
    print("The question is: " + question)
    print("\n--- Analysis ---")
    print(explanation)
    print("\n--- Conclusion ---")
    print(f"Character: {correct_character} => Final Answer Letter: {correct_option}")

solve_solaris_question()