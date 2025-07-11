def decode_word(encoded_word):
    """
    Decodes a word based on the discovered pattern.
    The pattern involves splitting the word into two, calculating a numeric value for each part,
    and then summing these values.
    """
    # 1. Split the word into two sub-strings based on character position (even/odd)
    s1 = encoded_word[0::2]
    s2 = encoded_word[1::2]

    # 2. Calculate the numeric value for each sub-string
    # The value is the sum of the letters' positions in the alphabet (a=0, b=1, ...)
    def calculate_value(s):
        total_value = 0
        for char in s:
            total_value += ord(char) - ord('a')
        return total_value

    val1 = calculate_value(s1)
    val2 = calculate_value(s2)

    # 3. The final decoded value is the sum of the two values
    final_value = val1 + val2

    # Print the equation as requested
    print(f"{val1} + {val2} = {final_value}")

    return final_value

# The last word to be decoded
last_word = "anhaeketn"

# Decode the word and print the equation
decode_word(last_word)