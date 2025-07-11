def decode_word_by_vowels(word):
    """
    Decodes a word by counting its vowels (a, e, i, o, u).
    The function prints the count of each vowel and the final sum as an equation.
    """
    vowels = "aeiou"
    vowel_counts = {v: 0 for v in vowels}
    
    for char in word.lower():
        if char in vowels:
            vowel_counts[char] += 1
            
    counts = []
    for v in vowels:
        counts.append(vowel_counts[v])
        
    total_vowels = sum(counts)
    
    # Build the equation string with each number
    equation_parts = [str(c) for c in counts]
    equation_str = " + ".join(equation_parts)
    
    print(f"Analyzing the word: '{word}'")
    print(f"Vowel counts (a, e, i, o, u): {counts}")
    print("The final equation is:")
    print(f"{equation_str} = {total_vowels}")

# The word to decode
target_word = "anhaeketn"
decode_word_by_vowels(target_word)