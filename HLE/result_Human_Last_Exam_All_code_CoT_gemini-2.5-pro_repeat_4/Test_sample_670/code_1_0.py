def generate_tate_expression():
  """
  Generates and prints the Tate-style program-expression for the given code.
  """
  # The variable 'a' is initialized to 1. This is the initial value for the recurrence.
  initial_value = 1

  # The update operation inside the loop is a = a * i.
  # This forms the recurrent part of the expression.
  # 'a' in the expression represents the value from the previous iteration (the bound variable).
  # 'i' is the loop counter.
  recurrent_expr = "a * i"

  # The μ (mu) operator is used to express the cycle (least fixed-point).
  # The PHI function takes the initial value and the recurrent value.
  # The bound variable, as requested, is 'a'.
  tate_expression = f"μ a . PHI({initial_value}, {recurrent_expr})"

  print(tate_expression)

generate_tate_expression()