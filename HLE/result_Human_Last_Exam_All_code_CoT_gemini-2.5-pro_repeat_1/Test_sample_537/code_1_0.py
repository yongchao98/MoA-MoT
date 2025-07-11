def decode_word(word):
  """
  This function decodes the word based on its length and prints the decoding equation.
  """
  length = len(word)
  
  # Create the equation string like "1 + 1 + ... + 1 = length"
  equation_parts = ["1"] * length
  equation_str = " + ".join(equation_parts)
  
  # Print the final equation with all its numbers
  print(f"The decoded value is calculated from the word's length:")
  print(f"'{word}' has {length} letters.")
  print(f"Final equation: {equation_str} = {length}")
  print(f"The decoded number is {length}.")

# The word to be decoded
final_word = "anhaeketn"
decode_word(final_word)