# The secret word and an equation can be derived from the image's right border.
# The border has a repeating sequence of 8 colored blocks, which can be mapped to numbers 0-7.
# Faint letters spelling the word "UZUME" are found within the blocks for numbers 1, 2, 3, 4, and 5.

# These are the numbers corresponding to the blocks containing the letters of the secret word.
num1 = 1  # 'U' in the blue block
num2 = 2  # 'Z' in the green block
num3 = 3  # 'U' in the cyan block
num4 = 4  # 'M' in the red block
num5 = 5  # 'E' in the magenta block

# The secret word is formed by the letters found.
secret_word = "UZUME"
print(f"The secret word found in the image is: {secret_word}")

# As requested, here is the equation formed by summing these numbers.
total = num1 + num2 + num3 + num4 + num5

# Print the final equation with each number.
print("The final equation is:")
print(f"{num1} + {num2} + {num3} + {num4} + {num5} = {total}")
