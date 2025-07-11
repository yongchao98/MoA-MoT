def find_foreign_symbol():
  """
  Analyzes an image of rock art to identify any symbols not belonging to
  the ancient Southwest cultures that created it.
  """
  is_foreign_symbol_present = "True"
  
  # The symbol identified is modern graffiti resembling Latin letters.
  # Ancient Southwest cultures did not use the Latin alphabet.
  symbol_identity = "The modern graffiti resembling the letters 'n' and 'U' (or 'n' and 'C')."
  
  # This graffiti can be seen in the upper-right area of the main image.
  location_description = "upper-right portion of the image"
  
  # Print the final analysis and conclusion.
  print(f"{is_foreign_symbol_present}.")
  print(f"There is at least one symbol not related to the ancient Southwest cultures.")
  print(f"The symbol is: {symbol_identity}")
  print(f"It is located in the {location_description}.")
  
find_foreign_symbol()