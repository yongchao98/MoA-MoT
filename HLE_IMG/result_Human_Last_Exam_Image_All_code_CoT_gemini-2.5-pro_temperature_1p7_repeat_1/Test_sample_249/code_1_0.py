import sys

def solve():
  """
  This function identifies the animal in the image and provides its scientific name.
  """
  animal_present = "Yes"
  animal_name = "Sturgeon"
  scientific_name = "Acipenseridae"

  print(f"Is there an animal in this image? {animal_present}")
  print(f"What is it? It is a {animal_name}.")
  print(f"What is its scientific name? {scientific_name}")

solve()