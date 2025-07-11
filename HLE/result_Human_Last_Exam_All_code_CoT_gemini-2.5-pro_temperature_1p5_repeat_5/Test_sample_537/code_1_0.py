def decode_word(word):
    """
    Decodes a word by splitting it in half, counting the vowels in each half,
    and finding the difference.
    """
    vowels = "aeiou"
    
    # Split the word into two halves
    split_index = len(word) // 2
    part1 = word[:split_index]
    part2 = word[split_index:]
    
    # Count vowels in the first half (v1)
    v1 = 0
    for char in part1:
        if char in vowels:
            v1 += 1
            
    # Count vowels in the second half (v2)
    v2 = 0
    for char in part2:
        if char in vowels:
            v2 += 1
            
    # Calculate the result
    result = v2 - v1
    
    # Print the equation as requested
    print(f"{v2} - {v1} = {result}")
    
    return result

# The last word to decode
last_word = "anhaeketn"

# Decode the word and print the final equation
final_answer = decode_word(last_word)

# The final answer is the result of the equation
# To follow the output format, this part is needed.
# print(f"<<<{final_answer}>>>")