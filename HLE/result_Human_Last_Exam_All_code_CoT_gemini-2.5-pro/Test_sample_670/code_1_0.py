def generate_tate_expression():
  """
  Generates and prints the Tate-style program-expression graph for the given loop.
  """

  # The recursive function is defined by the loop's behavior.
  # Let the state be (a, b), where 'a' maps to the program variable 'a'
  # and 'b' maps to the loop counter 'i'.
  
  # Recursive function body: "if b < n then f(a * b, b + 1) else a"
  # Lambda expression: "λ(a, b) . if b < n then f(a * b, b + 1) else a"
  
  # Using the mu-operator to define the recursive function f:
  # "μ f . λ(a, b) . if b < n then f(a * b, b + 1) else a"
  
  # The initial state is (a=1, i=1).
  # The final expression is the function applied to the initial state.
  
  expression = "(μ f . λ(a, b) . if b < n then f(a * b, b + 1) else a)(1, 1)"
  
  print(expression)

generate_tate_expression()