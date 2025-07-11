import math

def is_length_perfect_square(w):
  """
  This function decides the language L = {w | |w| is a perfect square}.
  It implements an algorithm that is guaranteed to halt for any input string w.
  """
  n = len(w)
  print(f"Input string: '{w}'")
  print(f"Length of the string is n = {n}.")
  print("Is n a perfect square? Let's check.")

  # Handle the case of an empty string
  if n == 0:
    print("Equation: 0 * 0 = 0")
    print(f"The length {n} is a perfect square. The string is in the language L.")
    return True
  
  i = 0
  while True:
    square = i * i
    print(f"Checking equation: {i} * {i} = {square}")
    
    if square == n:
      print(f"Success! The length {n} is a perfect square.")
      print("The string is in the language L.")
      return True
    
    if square > n:
      print(f"The square {square} is greater than the length {n}.")
      print(f"The length {n} is not a perfect square.")
      print("The string is NOT in the language L.")
      return False
      
    i += 1

# --- Main execution ---
# Let's test with a string of length 16, which is a perfect square.
input_string = "a" * 16

# Run the decider for the language L
is_length_perfect_square(input_string)

print("\n--------------------------\n")

# Let's test with a string of length 15, which is not a perfect square.
input_string_2 = "b" * 15
is_length_perfect_square(input_string_2)
