def reveal_hidden_word():
    """
    Solves a series of video game trivia questions to reveal a hidden word.
    The hidden word is formed by combining the second letter of each answer.
    """
    # List of answers to the trivia questions
    answers = [
        "OCCUPY",
        "KONAMI",
        "EARTHWORM JIM",
        "LEO",
        "MOUSE HAND"
    ]
    
    question_prompts = [
        "(1) Soviet Union / Moon action:",
        "(2) Received apology for anthem mix-up:",
        "(3) Game with 'Eat your vegetables' sound:",
        "(4) Pantera's companion (Latin origin):",
        "(5) Dennis Fong's 'MAGIC HAND' phrase:"
    ]

    hidden_word_letters = []
    
    print("Finding the hidden word by taking the second letter of each answer:")
    print("-" * 60)

    for i, answer in enumerate(answers):
        # The second letter is at index 1
        second_letter = answer[1]
        hidden_word_letters.append(second_letter)
        
        # Print the breakdown for each question
        print(f"Question {i+1}: {question_prompts[i]:<40} Answer: {answer:<15} -> Second Letter: {second_letter}")

    # Combine the letters to form the final word
    final_word = "".join(hidden_word_letters)
    
    print("-" * 60)
    print(f"The combined second letters form the hidden word: {final_word}")

# Run the function to reveal the word
reveal_hidden_word()