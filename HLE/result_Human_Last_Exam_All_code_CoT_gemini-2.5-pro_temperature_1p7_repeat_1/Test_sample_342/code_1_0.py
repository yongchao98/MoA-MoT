def calculate_Dg(g):
  """
  Calculates the g-th term of the sequence D_g.
  The formula D_g = 2**(g-1) * (2**g + 1) represents the number of
  even theta characteristics on a g-dimensional principally polarised
  abelian variety, which corresponds to the degree of the required cover.
  """
  # Check if g is a positive integer
  if not isinstance(g, int) or g < 1:
    raise ValueError("g must be a positive integer.")
  
  term = (2**(g - 1)) * (2**g + 1)
  return term

def generate_sequence():
  """
  Generates and prints the first 4 terms of the sequence D_g.
  """
  sequence = []
  for g in range(1, 5):
    dg_value = calculate_Dg(g)
    sequence.append(str(dg_value))
  
  # The problem asks to write the sequence. We will print the four numbers
  # as a comma-separated string. The prompt's instruction to "output
  # each number in the final equation" is interpreted as simply providing
  # the computed numbers of the sequence, as there is no literal equation here.
  print(", ".join(sequence))

generate_sequence()