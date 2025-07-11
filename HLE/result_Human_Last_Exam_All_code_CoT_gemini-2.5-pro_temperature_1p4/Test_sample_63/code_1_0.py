def solve_puzzle():
    # The man Kurt Vonnegut described as looking like a porcupine is the
    # German philosopher Arthur Schopenhauer. While his actual name doesn't
    # fit the constraint, a common title used for him does.

    first_word = "Pessimistic"
    second_word = "Philosopher"

    # The prompt requires an equation with numbers. We will create one
    # using the letter count of the two words in the answer.
    first_number = len(first_word)
    second_number = len(second_word)
    result = first_number + second_number

    # Print the name, which is the primary answer to the user's question.
    print(f"The name for the man described is: {first_word} {second_word}")

    # Print the full equation, which includes each number as requested.
    print(f"A novelty equation using the letter count is: {first_number} + {second_number} = {result}")

solve_puzzle()