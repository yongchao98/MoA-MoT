def sum_digits(n):
  """Calculates the sum of the digits of a number."""
  s = 0
  for digit in str(n):
    s += int(digit)
  return s

def find_next_number_in_sequence(start_num):
  """
  Finds the next number in the sequence starting from start_num.
  The rule is that the number must be divisible by 6 and also be a
  Harshad number (divisible by the sum of its digits).
  """
  num = start_num
  while True:
    # We only need to check numbers greater than start_num.
    # We can increment by 6 to ensure the number is always divisible by 6.
    num += 6
    
    s = sum_digits(num)
    
    # Check if it's a Harshad number
    if num % s == 0:
      return num

# The last number in the given sequence
last_number = 2352

# Find the next number
next_number = find_next_number_in_sequence(last_number)

# The "equation" demonstrates that the found number satisfies the Harshad property.
# The divisibility by 6 is guaranteed by the search method.
print(f"The last number is {last_number}.")
print(f"The next number in the sequence is {next_number}.")
print("The rule is that the number must be divisible by the sum of its digits (a Harshad number) and also be divisible by 6.")
print(f"For {next_number}, the sum of its digits is {sum_digits(next_number)}.")
print(f"The final equation based on the Harshad property is:")
print(f"{next_number} % ({' + '.join(str(next_number))}) = {next_number % sum_digits(next_number)}")