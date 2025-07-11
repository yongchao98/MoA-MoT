def find_substituents():
  """
  Identifies the substituents at the numbered positions in the product molecule
  based on standard chemical structure notation.
  """
  substituents = {
      1: "CH3",
      2: "CH3",
      3: "H",
      4: "CH3",
      5: "H"
  }

  # Print the result in the specified format
  # We use a list to ensure the output order is as requested.
  output_list = []
  for i in range(1, 6):
      output_list.append(f"{i} = {substituents[i]}")

  print(", ".join(output_list))

find_substituents()