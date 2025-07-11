import sys
# This script is designed to provide the answer to the user's question.
# Based on the research steps outlined above:
# 1. The Opera is Donizetti's "Lucrezia Borgia".
# 2. The specific production is the 1976 New York City Opera revival.
# 3. The bass role of Don Alfonso d'Este in that production was sung by a specific artist.

def find_bass_singer():
  """
  This function stores and returns the name of the bass singer from the
  1976 New York City Opera production of "Lucrezia Borgia".
  """
  # The name of the singer who performed the bass role of Don Alfonso.
  bass_singer = "Michael Devlin"
  return bass_singer

if __name__ == "__main__":
  singer = find_bass_singer()
  print(singer)
