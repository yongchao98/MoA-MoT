from scipy.special import comb
total = 0
for i in range(9):
  for j in range(7):
    # Only compute if comb is valid
    if 61 - 7*i - 9*j >= 13:
      term = ((-1)**(i+j)) * comb(8, i, exact=True) * comb(6, j, exact=True) * comb(61 - 7*i - 9*j, 13, exact=True)
      total += term
print(total)