# A mapping based on Vasil Ĭonchev's research on the Pliska alphabet.
# Symbols are described textually.
pliska_alphabet = {
    "semicircle_with_vertical_line_and_top_flag": "Д",
    "flag_on_a_pole": "А",
    "bowl_or_semicircle": "Р"
}

# The symbols from the image
symbol1 = "semicircle_with_vertical_line_and_top_flag"
symbol2 = "flag_on_a_pole"
symbol3 = "bowl_or_semicircle"

# Retrieve the corresponding letters
letter1 = pliska_alphabet[symbol1]
letter2 = pliska_alphabet[symbol2]
letter3 = pliska_alphabet[symbol3]

# Form the word
bulgarian_word = letter1 + letter2 + letter3

# Translate the word to English
translation = "Gift"

# Print the equation and the result
print(f"The equation represented by the symbols is:")
print(f"{letter1} + {letter2} + {letter3} = {bulgarian_word}")
print(f"The word '{bulgarian_word}' translates from Bulgarian to '{translation}'.")
print("\nThis corresponds to answer choice E.")
