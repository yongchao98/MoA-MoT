def decode_word_to_equation(word):
    """
    Decodes a word by converting each letter to its alphabetical position
    (A=1, B=2, ...), and then prints the summation equation.
    """
    word = word.lower()
    numbers = []
    total_sum = 0
    
    for char in word:
        if 'a' <= char <= 'z':
            # Calculate the value (a=1, b=2, ...)
            value = ord(char) - ord('a') + 1
            numbers.append(value)
            total_sum += value
            
    # Build the equation string
    equation_parts = [str(n) for n in numbers]
    equation_str = " + ".join(equation_parts)
    
    # Print the full equation including the result
    print(f"{equation_str} = {total_sum}")
    
    # Return the final sum as the decoded value
    return total_sum

# The last word to be decoded
last_word = "anhaeketn"

# Decode the word and get the final answer
final_answer = decode_word_to_equation(last_word)

# The final answer is the sum from the equation
# The print statement above will show the step-by-step equation
# Now, we will output the final answer in the required format.
# print(f"<<<{final_answer}>>>")