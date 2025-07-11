import sys

def solve():
  """
  This function provides the formula for P(n) based on the asymptotic analysis.
  The derivation is outlined in the thinking steps. The core idea is to find
  the asymptotic series for Q(n) and determine which terms P(n) must represent
  to achieve the desired error bound.
  """
  # The problem asks for the formula for P(n), where L = ln(n).
  # The derivation shows that P(n) must cancel the O((L/n)^2) and O((L/n)^3) terms
  # in the asymptotic expansion of Q(n) / (A * n^(L/2)).
  # These terms are found to be:
  # C2/n^2 = (3*L^2 + 2*L - 2) / (24*n^2)
  # C3/n^3 = (L^3 + 2*L^2 - 2*L) / (48*n^3)
  # So, P(n) = C2/n^2 + C3/n^3.

  term1_num_L2_coeff = 3
  term1_num_L_coeff = 2
  term1_num_const = -2
  term1_den = 24

  term2_num_L3_coeff = 1
  term2_num_L2_coeff = 2
  term2_num_L_coeff = -2
  term2_den = 48

  # Note: The output format uses standard mathematical notation, not Python syntax.
  print(f"P(n) = ({term1_num_L2_coeff}*L^2 + {term1_num_L_coeff}*L + {term1_num_const})/({term1_den}*n^2) + ({term2_num_L3_coeff}*L^3 + {term2_num_L2_coeff}*L^2 + {term2_num_L_coeff}*L)/({term2_den}*n^3)")

solve()