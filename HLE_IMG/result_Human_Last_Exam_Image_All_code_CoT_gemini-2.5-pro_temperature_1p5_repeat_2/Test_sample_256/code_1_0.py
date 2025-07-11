def solve_image_puzzle():
  """
  This function analyzes the rock art image and provides the answer.
  """
  is_foreign_symbol_present = True
  foreign_symbol_description = "the number 70"

  # The main body of art is from ancient Southwest cultures.
  # However, in the upper right section, modern graffiti in the form of Arabic numerals is visible.
  # These numerals are not part of the original culture's symbology.
  
  answer = f"{is_foreign_symbol_present}.\nOn the upper right side of the image, it is possible to distinguish two symbols that form {foreign_symbol_description}. These Arabic numerals are not part of the symbology of the ancient Southwest cultures that created the rest of the pictographs. This is likely a later addition or graffiti."

  print(answer)

solve_image_puzzle()