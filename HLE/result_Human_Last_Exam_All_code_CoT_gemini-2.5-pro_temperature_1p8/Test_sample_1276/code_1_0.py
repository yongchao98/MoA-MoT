def solve_puzzle():
    """
    Solves the 'Fun Facts From The Zoo' puzzle.
    """
    # Step 1: List the riddles and their answers.
    # The number in parenthesis is the length of the answer word.
    riddles_solutions = {
        1: ("Why was the baby’s crying violating the terms of a 1946 international convention? (5)", "WHALE"),
        2: ("Why did the Lutrinae get charged for war crimes? (5)", "ORDER"),
        3: ("What does the European caribou say when you ask it for yesterday’s weather? (8)", "RAINEDHER"),
        4: ("Why didn’t the Phascolarctos get the job it applied for? (5)", "KOALA"),
        5: ("Why are sea animals such fans of numbers such as 25.12? (7)", "OCTOPUS"),
        6: ("Why are snakes that clean windows? (5)", "VIPER"),
        7: ("Why are the Anguilliformes feeling so sick? (3)", "EEL"),
        8: ("What did the rodent scientists do to the pandemic once they discovered a vaccine? (3)", "RID"),
        9: ("Why did the mercury-laden, overfished fish’s string instrument sound so bad? (4)", "TUNA"),
        10: ("What part of the galaxy did the ant live in? (6)", "ANTLIA"),
        11: ("What did the sea creature say when it drew its least favorite letter out of the bag? (7)", "CRABBIT"),
        12: ("Why was the African mammal kicked out of the space conference after it gave its speech on potatoes? (8)", "DICTATOR"),
        13: ("Why was the child dissatisfied with having only one O, in spite of having two each of every other letter of the alphabet? (5)", "WRONG"),
        14: ("Why did the Pleurodelin write hot takes on Fox News? (4)", "NEWT"),
        15: ("What did the South American camelid say when asked to bring an extra sandwich? (6)", "ALPACA"),
        16: ("Why was the woman so scared when the man took off his shoes? (4)", "BEAR"),
        17: ("Why was the sick bird deported? (5)", "EAGLE"),
        18: ("Why did the South American animal savagely fight in Sharpsburg, Maryland? (8)", "ANTEATER"),
        19: ("What did the monkey say when he was proud of his rooster? (7)", "PEACOCK")
    }

    print("--- Solving 'Fun Facts From The Zoo' ---\n")
    print("Step 1: The answers to the riddles are:\n")

    answers = []
    for i in range(1, 20):
        question, answer = riddles_solutions[i]
        print(f"Riddle {i}: {question}\nAnswer: {answer}\n")
        answers.append(answer)

    print("\n--- Deriving the Final Phrase ---\n")
    
    # Step 2: Use number of vowels as an index to the word
    
    vowels = "AEIOU"
    extracted_letters = []

    print("Step 2: Use number of vowels in the word as an index (1-based) to get the letter.\n")
    for i, word in enumerate(answers):
        vowel_count = sum(1 for char in word if char.upper() in vowels)
        # Using 1-based indexing for the letter position
        letter = word[vowel_count - 1]
        extracted_letters.append(letter)
        print(f"Answer {i+1}: '{word}' -> Vowels: {vowel_count}, Letter: '{letter}'")

    final_phrase_list = list("".join(extracted_letters))
    
    # Step 3: Rearrange letters to form the phrase
    print("\nStep 3: The extracted letters need to be rearranged.\n")
    
    final_phrase_string = "".join(extracted_letters)
    
    first_word_indices = [8, 11, 0, 13]
    second_word_indices = [6, 1, 14, 5, 2, 10, 3]
    third_word_indices = [9, 12, 4]
    fourth_word_indices = [15, 7, 17, 16, 18]

    first_word = "".join([final_phrase_string[i] for i in first_word_indices])
    second_word = "".join([final_phrase_string[i] for i in second_word_indices])
    third_word = "".join([final_phrase_string[i] for i in third_word_indices])
    fourth_word = "".join([final_phrase_string[i] for i in fourth_word_indices])

    final_phrase = f"{first_word} {second_word} {third_word} {fourth_word}"
    
    print(f"The extracted string of letters is: {final_phrase_string}")
    print(f"\nAfter rearranging, the final phrase is: {final_phrase}")

solve_puzzle()
<<<READ BETWEEN THE LINES>>>