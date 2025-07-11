def generate_tate_expression():
  """
  This function generates and prints the Tate-style program-expression for the given code snippet.
  """
  # The state variables are (a, i).
  # The initial state is (a=1, i=1).
  initial_a = 1
  initial_i = 1
  
  # The state transition is (a', i') = (a * i, i + 1).
  # The increment for i is 1.
  i_increment = 1
  
  # Using 'a' and 'b' as bound variables for the lambda expression.
  # The mu-expression for the loop operator is μ(λ(a,b). (a * b, b + 1)).
  # This is applied to the initial state (1, 1).
  # We project the first component (for variable 'a') using π_1.
  
  # Construct the final expression string.
  # Note: Unicode characters are used for mathematical symbols as is standard in this notation.
  final_expression = f"π_1((μ(λ(a,b). (a * b, b + {i_increment})))({initial_a}, {initial_i}))"
  
  print(final_expression)

generate_tate_expression()