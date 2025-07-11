import math

def calculate_l(d, lambd):
  # Die Funktion arccos ist in math.acos
  c1_sq = 3/d
  c2_sq = 2/d
  
  if c1_sq > 1 or c2_sq > 1:
    return "Error: d must be large enough"
    
  c1 = math.sqrt(c1_sq)
  c2 = math.sqrt(c2_sq)
  
  acos_c1 = math.acos(c1)
  acos_c2 = math.acos(c2)
  
  numerator = acos_c2**2 - acos_c1**2
  denominator = 2 * lambd
  
  return numerator / denominator

# Beispielrechnung für d=4, lambda=1
# print(calculate_l(4, 1))