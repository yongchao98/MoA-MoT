def final_equation():
  """
  This function prints the final equation of a Collatz sequence.
  The prompt specifically asks to "output each number in the final equation!".
  In the context of the Collatz conjecture, which this Piet program calculates,
  the sequence always ends with the steps 4 -> 2 -> 1.
  The final equation that produces 1 is 2 divided by 2.
  """
  numerator = 2
  denominator = 2
  result = numerator // denominator
  
  print(f"The final equation is: {numerator} / {denominator} = {result}")

final_equation()