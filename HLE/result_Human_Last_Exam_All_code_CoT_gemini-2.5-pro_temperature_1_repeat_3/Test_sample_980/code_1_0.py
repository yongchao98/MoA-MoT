import functools
import operator
import re

def calculate_word_value(word):
    """Calculates the value of a word by XORing the ASCII values of its characters."""
    if not word:
        return 0
    # The 'equation' is the XOR sum of the character ordinals.
    char_ords = [ord(c) for c in word]
    value = functools.reduce(operator.xor, char_ords)
    return value

def solve_sequence():
    """
    Solves the sequence puzzle by finding the generating rule and applying it.
    """
    # The source text is the full prompt itself.
    source_text = """You are a helpful AI assistant. Solve tasks using your coding skills.
    In the following cases, suggest python code (in a python coding block) or shell script (in a sh coding block) for the user to execute.
    Don't include multiple code blocks in one response, only include one in the response. Do not ask users to copy and paste the result. Instead, use 'print' function for the output when relevant.
    Think the task step by step if you need to. If a plan is not provided, explain your plan first. You can first output your thinking steps with texts and then the final python code.
    Remember in the final code you still need to output each number in the final equation!
    Start the python block with ```python

    What is the next number in the sequence:"""

    # Tokenize the text into words, keeping only alphanumeric characters.
    words = re.findall(r'\b\w+\b', source_text)

    given_sequence = [111, 142, 111, 41, 67, 67, 67, 93, 111, 111, 62, 62, 111, 111, 36, 36, 49, 155, 49, 62, 49, 49, 62, 62, 10, 36, 36, 36, 124, 124, 124, 36, 124]
    
    # Find the starting index of the sequence in the text.
    start_index = -1
    for i in range(len(words)):
        # The word "Don't" is the start of the sequence.
        if words[i] == "Dont":
            start_index = i
            break
            
    if start_index == -1:
        print("Error: Could not find the start of the sequence in the source text.")
        return

    # The word that generates the next number is the one after the end of the given sequence.
    next_word_index = start_index + len(given_sequence)
    next_word = words[next_word_index]

    # Calculate the next number
    next_number = calculate_word_value(next_word)

    # Per instructions, output the final equation for the next number.
    equation_parts = [f"ord('{c}')" for c in next_word]
    equation_str = " ^ ".join(equation_parts)
    
    ord_values = [str(ord(c)) for c in next_word]
    ord_equation_str = " ^ ".join(ord_values)

    print(f"The next word in the sequence is '{next_word}'.")
    print(f"The equation to find the next number is: {equation_str}")
    print(f"Which evaluates to: {ord_equation_str} = {next_number}")
    print(f"\nThe next number in the sequence is: {next_number}")


solve_sequence()