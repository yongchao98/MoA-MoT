import unicodedata

def get_radical_name_in_hepburn(character):
  """
  This function provides the Hepburn romanization of the radical for a specific kanji.
  This is a simplified approach using a predefined dictionary for the requested character.
  """
  kanji_radicals = {
      '謄': 'Nikuzuki', # The radical is 肉 (niku), which appears as 月 (tsuki) on the left. This form is called nikuzuki (にくづき).
      '媛': 'Onnahen'   # Example from the prompt.
  }
  
  radical_name = kanji_radicals.get(character, "Radical not found in this simple dictionary.")
  print(radical_name)

# The character in question is '謄'.
character_to_find = '謄'
get_radical_name_in_hepburn(character_to_find)