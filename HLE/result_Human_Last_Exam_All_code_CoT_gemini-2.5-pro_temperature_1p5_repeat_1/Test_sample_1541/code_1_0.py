import sys

def solve_xenia_sum():
  """
  This function reads a string representing a sum of numbers 1, 2, 3, and 4,
  and prints the sum with the numbers sorted in non-decreasing order.
  """
  
  # Read the entire line from standard input and remove any trailing whitespace.
  input_sum = sys.stdin.readline().strip()
  
  # Split the string by the '+' character to get a list of the numbers (as strings).
  # For example, "1+3+2+1" becomes ['1', '3', '2', '1'].
  summands = input_sum.split('+')
  
  # Sort the list of numbers. Python's default sort for strings works
  # correctly here because all numbers are single digits.
  # ['1', '3', '2', '1'] becomes ['1', '1', '2', '3'].
  summands.sort()
  
  # Join the sorted list of numbers back into a single string,
  # with '+' as the separator between them.
  # ['1', '1', '2', '3'] becomes "1+1+2+3".
  # This final string represents the full equation with each number included.
  sorted_sum = "+".join(summands)
  
  # Print the final result.
  print(sorted_sum)

# Execute the function to solve the problem.
solve_xenia_sum()