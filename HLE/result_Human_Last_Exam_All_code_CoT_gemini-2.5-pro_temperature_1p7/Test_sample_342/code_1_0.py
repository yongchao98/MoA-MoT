def calculate_dg(g):
  """
  Calculates the g-th term of the sequence D_g.
  D_g = 2**(g-1) * (2**g + 1)
  """
  return 2**(g-1) * (2**g + 1)

# Calculate the first 4 terms of the sequence, starting from g=1.
g_values = [1, 2, 3, 4]
sequence = []

print("The formula for D_g is: D_g = 2**(g-1) * (2**g + 1)")
print("-" * 30)

for g in g_values:
  dg_value = calculate_dg(g)
  sequence.append(dg_value)
  print(f"For g = {g}, D_{g} = 2**({g-1}) * (2**{g} + 1) = {2**(g-1)} * {(2**g + 1)} = {dg_value}")

print("-" * 30)
print("The first 4 terms of the sequence D_g are:", ", ".join(map(str, sequence)))