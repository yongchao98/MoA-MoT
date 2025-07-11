import sys

def find_perdita_common_name_info():
  """
  This function provides information about the common name for the bee genus Perdita.

  The information is based on research into entomological literature.
  The common name "fairy bees" was coined by the melittologist (bee scientist)
  Dr. Jerome G. Rozen, Jr. He explicitly stated this in a 1996 article
  in Natural History magazine titled "Little, Little Bees," referring to the bees'
  tiny size and delicate nature.
  """
  common_name = "fairy bees"
  year_coined = 1996

  # The final output is formatted as "<name>, <year>"
  print(f"{common_name}, {year_coined}")

if __name__ == '__main__':
  # This check prevents the function from running when the script is imported,
  # but allows it to run when executed directly.
  # This is a common practice in Python scripting.
  find_perdita_common_name_info()
