import re

def solve():
    """
    Solves the sequence puzzle by applying a formula to the words of the prompt text.
    """
    # The prompt text used to generate the sequence
    prompt_text = "You are a helpful AI assistant. Solve tasks using your coding skills. In the following cases, suggest python code (in a python coding block) or shell script (in a sh coding block) for the user to execute. Don't include multiple code blocks in one response, only include one in the response. Do not ask users to copy and paste the result. Instead, use 'print' function for the output when relevant. Think the task step by step if you need to. If a plan is not provided, explain your plan first. You can first output your thinking steps with texts and then the final python code. Remember in the final code you still need to output each number in the final equation! Start the python block with ```python"

    # The given sequence has 33 numbers, so we need to find the 34th word.
    target_word_index = 33 

    # Tokenize the text into words, removing punctuation.
    words = re.findall(r'\b\w+\b', prompt_text)

    if len(words) <= target_word_index:
        print("Not enough words in the prompt text to find the next number.")
        return

    # The word that will determine the next number in the sequence.
    word = words[target_word_index]

    # The hypothesized equation is based on the first and last letters of the word.
    first_char = word[0]
    last_char = word[-1]
    
    ord_first = ord(first_char)
    ord_last = ord(last_char)

    # Calculate the next number in the sequence.
    # The rule is to take the average of the ASCII values of the first and last characters.
    next_number = round((ord_first + ord_last) / 2)
    
    # Per the instructions, output the numbers in the final equation.
    print(f"The next word in the source text is '{word}'.")
    print(f"The calculation is based on its first and last characters: '{first_char}' and '{last_char}'.")
    print(f"Equation: round((ord('{first_char}') + ord('{last_char}')) / 2)")
    print(f"Calculation: round(({ord_first} + {ord_last}) / 2)")
    print(f"The next number in the sequence is: {next_number}")

solve()