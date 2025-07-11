def solve():
  """
  This function identifies and prints the letters corresponding to the statements
  that describe necessarily Lebesgue integrable functions.
  """
  # Based on mathematical analysis:
  # F: A continuous function on [a,b] is bounded and measurable, so integrable on a set of finite measure.
  # H: If a function's absolute value is integrable, it's measurable and hence satisfies the definition of Lebesgue integrability.
  # I: Same as H, but on a specific finite measure space [a,b].
  # L: This is the explicit definition of a Lebesgue integrable function.
  # M: A bounded continuous function is measurable, and its integral over a finite interval (a,b) is finite.
  
  correct_choices = "FHILM"
  print(correct_choices)

solve()