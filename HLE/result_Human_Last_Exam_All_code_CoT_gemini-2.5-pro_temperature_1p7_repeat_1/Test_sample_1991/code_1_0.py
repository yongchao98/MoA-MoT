def generate_descriptive_phrase():
  """
  This function generates and prints a fourteen-syllable phrase.
  The phrase describes abnormal cell death causing cancer, with rhyme
  and scansion matching a famous Disney song title.
  """
  # The phrase is broken down for syllable analysis:
  # "A-pop-to-sis" (4) + "Reg-u-la-tion" (4) = 8 syllables
  # "Which" (1) + "Is" (1) + "Quite" (1) + "A-tro-cious" (3) = 6 syllables
  # Total = 14 syllables
  
  part_one = "Apoptosis Regulation"
  part_two = "Which Is Quite Atrocious"
  
  final_phrase = f"{part_one} {part_two}"
  
  print(final_phrase)

generate_descriptive_phrase()