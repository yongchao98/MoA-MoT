def solve():
  """
  Determines for six classes of preordered sets whether they always have
  a maximal element (Y), never do (N), or if it depends (D).
  """

  # Case A: H-free graphs. Depends on H (e.g., K1 vs K2).
  answer_A = 'D'

  # Case B: Finite discrete subset of R. Assuming non-empty, always has a maximum.
  answer_B = 'Y'

  # Case C: Countable discrete subset of R. Depends on the set (e.g., N vs {-n}).
  answer_C = 'D'

  # Case D: Uncountable discrete subset of R. This class is empty, so the
  # statement "all have a maximal element" is vacuously true.
  answer_D = 'Y'

  # Case E: Sequences with subsequence order. A universal sequence exists and is a greatest/maximal element.
  answer_E = 'Y'

  # Case F: Sequences with supersequence order. Constant sequences are maximal elements.
  answer_F = 'Y'

  final_answer = answer_A + answer_B + answer_C + answer_D + answer_E + answer_F
  print(final_answer)

solve()