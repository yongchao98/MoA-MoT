def calculate_Dg(g):
  """
  Calculates the number of odd theta characteristics for a PPAV of dimension g.
  This corresponds to the degree D_g.
  The formula is D_g = 2^(g-1) * (2^g - 1).
  """
  return (2**(g - 1)) * (2**g - 1)

# Calculate the first 4 terms of the sequence D_g, starting from g=1
g_values = [1, 2, 3, 4]
Dg_sequence = [calculate_Dg(g) for g in g_values]

# Print the sequence
print("The sequence of the first 4 terms of D_g is:")
# The problem asks to output each number in the final equation.
# Let's format it as a sequence.
print(', '.join(map(str, Dg_sequence)))
