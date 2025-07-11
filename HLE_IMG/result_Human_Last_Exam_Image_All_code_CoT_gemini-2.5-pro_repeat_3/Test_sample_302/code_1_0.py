import sys
import io

# Set stdout to utf-8 to ensure Chinese characters are printed correctly.
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

def get_seal_characters():
  """
  Identifies and returns the characters on the seal in simplified Chinese.
  The reading order is top-to-bottom, right-to-left.
  """
  # Characters identified from the seal image.
  top_right = "景"
  bottom_right_simplified = "园" # Traditional is 園
  top_left = "不"
  bottom_left_simplified = "识" # Traditional is 識

  # Combine the characters in the correct reading order.
  seal_text = top_right + bottom_right_simplified + top_left + bottom_left_simplified
  return seal_text

# Print the identified characters.
print(get_seal_characters())