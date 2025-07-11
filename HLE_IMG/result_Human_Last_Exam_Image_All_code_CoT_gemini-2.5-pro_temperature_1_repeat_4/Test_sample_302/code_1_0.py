def get_seal_characters():
  """
  This function returns the simplified Chinese characters found on the seal.
  The characters are read from top to bottom, right to left.
  """
  # The characters identified on the seal in traditional Chinese are 書, 畫, 不, 如.
  # Converted to simplified Chinese, they are:
  char1 = "书"  # shū (calligraphy)
  char2 = "画"  # huà (painting)
  char3 = "不"  # bù (not)
  char4 = "如"  # rú (as good as)
  
  # The phrase reads "shū huà bù rú"
  full_phrase = f"{char1}{char2}{char3}{char4}"
  
  print(f"The simplified Chinese characters on the seal are: {full_phrase}")

get_seal_characters()