def solve():
  """
  This function identifies and prints the simplified Chinese characters on the seal.
  The characters are read from top to bottom, right to left.
  """
  # The characters in traditional script are 壽康年豐.
  # Converted to simplified Chinese, they are 寿康年丰.
  character_1 = "寿"  # simplified of 壽 (shòu - longevity)
  character_2 = "康"  # 康 (kāng - health)
  character_3 = "年"  # 年 (nián - year)
  character_4 = "丰"  # simplified of 豐 (fēng - abundant)

  # The phrase is read right column (top to bottom), then left column (top to bottom).
  final_phrase = f"{character_1}{character_2}{character_3}{character_4}"

  print(f"The simplified Chinese characters on the seal read: {final_phrase}")
  # The final output is an auspicious phrase meaning "longevity, health, and a year of abundance".

solve()