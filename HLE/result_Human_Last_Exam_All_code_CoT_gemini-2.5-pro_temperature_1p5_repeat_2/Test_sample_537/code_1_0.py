def decode_word(word):
    """
    Decodes a word into a number by counting its vowels and consonants
    and returning their product.
    """
    vowels = "aeiou"
    vowel_count = 0
    consonant_count = 0
    for char in word.lower():
        # Ensure we only count alphabetic characters
        if 'a' <= char <= 'z':
            if char in vowels:
                vowel_count += 1
            else:
                consonant_count += 1
    return vowel_count, consonant_count

def solve_puzzle():
    """
    Solves the word puzzle by decoding each word and finding the final answer.
    """
    words = ["utpshtheas", "fkuhu", "numhkatnatu", "anhaeketn"]
    decoded_numbers = []

    for word in words:
        v, c = decode_word(word)
        decoded_numbers.append(v * c)

    n1, n2, n3, n4 = decoded_numbers

    # The prompt requires printing the equation formed by the first three numbers.
    print(f"{n1} + {n2} = {n3}")
    
    # The final answer is the decoding of the last word.
    # We will print it in the required format. This part will not be visible
    # in the final output block, but it's here for logical completeness.
    # print(f"<<<{n4}>>>")


solve_puzzle()
# The value for the last word 'anhaeketn' is 4 vowels * 5 consonants = 20.
# The code above will print the equation. The final answer is 20.