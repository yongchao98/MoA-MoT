import sys

# This script generates a fourteen-syllable phrase describing abnormal cell death
# in the Luys' body (subthalamic nucleus) that causes cancer.
# The phrase's scansion and rhyme match the Disney song title "Supercalifragilisticexpialidocious".
# Each word is capitalized as requested.

def generate_phrase():
  """
  Constructs and returns the final phrase.
  """
  # The phrase is composed of medical and descriptive terms to fit the 14-syllable count and theme.
  # "Subthalamic" refers to the Body of Luys.
  # "Cellular Necrosis" refers to abnormal cell death.
  # "Is Atrocious" implies the cancerous and harmful nature of the condition.
  phrase = "Subthalamic Cellular Necrosis Is Atrocious"
  return phrase

# Print the final phrase to the console.
final_phrase = generate_phrase()
print(final_phrase)

# For clarity, let's also output the source of inspiration.
# This part is for context and not part of the primary answer.
# song_title = "Supercalifragilisticexpialidocious"
# print(f"\nInspired by the scansion of: {song_title}", file=sys.stderr)
