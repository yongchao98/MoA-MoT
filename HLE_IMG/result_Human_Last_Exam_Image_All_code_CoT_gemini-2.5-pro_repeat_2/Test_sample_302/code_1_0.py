import sys
import io

# Set stdout to utf-8 to handle Chinese characters
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

def get_seal_characters():
  """
  This function identifies and returns the characters on the seal.
  The reading order is top-to-bottom, right-to-left.
  """
  # Traditional characters identified from the seal
  char1_trad = "長"
  char2_trad = "風"
  char3_trad = "不"
  char4_trad = "驚"

  # Simplified Chinese equivalents
  char1_simp = "长"
  char2_simp = "风"
  char3_simp = "不"
  char4_simp = "惊"

  return [char1_simp, char2_simp, char3_simp, char4_simp]

def main():
  """
  Main function to print the seal content.
  """
  characters = get_seal_characters()
  phrase = "".join(characters)
  
  print("The characters on the seal are read from top-to-bottom, right-to-left.")
  print(f"The phrase in simplified Chinese is: {phrase}")
  print("The individual characters are:")
  print(f"1. (Top-Right): {characters[0]}")
  print(f"2. (Bottom-Right): {characters[1]}")
  print(f"3. (Top-Left): {characters[2]}")
  print(f"4. (Bottom-Left): {characters[3]}")

if __name__ == "__main__":
  main()