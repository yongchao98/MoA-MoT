import math

def number_of_cool_strings(n):
  """
  Calculates the number of cool strings of maximal length for n symbols.
  The formula is (n!)^2.
  """
  if n < 0:
    return 0
  
  factorial_n = math.factorial(n)
  return factorial_n * factorial_n

def main():
  """
  Prints the number of cool strings for n=1 to 5.
  """
  print("The number of cool strings of maximal length for n symbols appears to follow the formula (n!)^2.")
  print("Here are the calculated values for n from 1 to 5:")
  for i in range(1, 6):
    n_factorial = math.factorial(i)
    result = n_factorial**2
    print(f"For n={i}: C({i}) = ({i}!)**2 = ({n_factorial})**2 = {result}")

if __name__ == "__main__":
  main()