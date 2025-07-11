import sys

def find_common_name_info():
  """
  This function provides the common name and the year it was noted
  in a scientific publication for the bee genus Perdita.

  The genus Perdita contains some of the world's smallest bees, and due to their
  tiny, delicate nature, they are often referred to by the common name "Fairy bees".

  While the exact origin of this informal name is difficult to pinpoint, one of the
  earliest documented mentions in a scientific paper appears in a 1996 publication
  by entomologist Bryan N. Danforth. In his paper "Nesting biology of the bee Perdita
  (Alloperdita) semicaerulea", he notes the "...small size of these bees (often
  referred to as 'fairy bees')...", indicating the name was already in use by that time.
  """
  common_name = "Fairy bees"
  year_documented = 1996

  # The final output needs to be in the format "<name>, <year>"
  # The print statement below constructs this string.
  print(f"{common_name}, {year_documented}")

if __name__ == "__main__":
  find_common_name_info()