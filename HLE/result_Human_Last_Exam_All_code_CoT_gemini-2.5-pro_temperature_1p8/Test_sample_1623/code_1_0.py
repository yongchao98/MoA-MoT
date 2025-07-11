import math

def sgn(x):
  """
  Helper function to compute the sign of a number.
  """
  if x > 0:
    return 1
  elif x < 0:
    return -1
  else:
    return 0

# Step 1: Define grid parameters from the problem description.
# The grid size is 5x5.
n = 5

# The circle 'o' positions are (1,1), (2,2), (3,3), (4,4), (5,5).
# In 0-based indexing, this corresponds to the permutation pi_O where the
# circle in column i is in row i.
pi_O = [0, 1, 2, 3, 4]

# The cross 'x' positions are (1,4), (2,5), (3,1), (4,2), (5,3).
# In 0-based indexing, this gives the permutation pi_X:
# Column 0 (was 1) has a cross in row 3 (was 4).
# Column 1 (was 2) has a cross in row 4 (was 5).
# Column 2 (was 3) has a cross in row 0 (was 1).
# Column 3 (was 4) has a cross in row 1 (was 2).
# Column 4 (was 5) has a cross in row 2 (was 3).
pi_X = [3, 4, 0, 1, 2]

# Step 2: Calculate the writhe w(G).
writhe = 0
for i in range(n):
  for j in range(i + 1, n):
    # Calculate the term for the pair (i,j)
    term = (pi_O[i] - pi_O[j]) * (pi_X[i] - pi_X[j])
    # Add its sign to the writhe
    writhe += sgn(term)

# Step 3: Calculate the Thurston-Bennequin number.
tb_number = writhe - n

# Step 4: Print the final calculation.
print("The Thurston-Bennequin number is calculated as w - n.")
print(f"Calculated writhe w = {writhe}")
print(f"Grid size n = {n}")
print(f"The final equation is: {writhe} - {n} = {tb_number}")
print(f"The maximal Thurston-Bennequin number is {tb_number}.")
