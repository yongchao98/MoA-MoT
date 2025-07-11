def find_next_number(n):
  """
  Calculates the next number in the sequence by summing the squares of the digits.
  This function also prints the equation used for the calculation as requested.
  """
  s_num = str(n)
  digits = [int(d) for d in s_num]
  
  # Build the equation string to display it
  equation_parts = []
  for d in digits:
    equation_parts.append(f"{d}^2")
  
  equation_str = " + ".join(equation_parts)
  
  # Calculate the result
  result = sum(d**2 for d in digits)
  
  # Print the equation, showing how each digit is used
  print(f"The calculation for the next number is based on the last number, {n}.")
  print(f"The equation is: {equation_str} = {result}")

# The last number in the given sequence
last_number = 124

find_next_number(last_number)