def generate_equation():
  """
  This function generates and prints an equation based on the answer to the riddle.
  The answer is "Felt hat". The numbers are derived from the alphabetical
  position of the letters in the words.
  """
  # Values are pre-calculated based on the plan:
  # Felt: F(6) + E(5) + L(12) + T(20) = 43
  # hat: H(8) + A(1) + T(20) = 29
  num1 = 43
  num2 = 29
  result = num1 + num2

  # Print the final equation, showing each number.
  print(f"{num1} + {num2} = {result}")

generate_equation()