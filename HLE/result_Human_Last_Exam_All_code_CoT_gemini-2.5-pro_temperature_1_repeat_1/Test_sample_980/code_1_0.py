def sum_digits(n):
  """Calculates the sum of the digits of a number."""
  s = 0
  for digit in str(n):
    s += int(digit)
  return s

# The last triplet in the sequence is (124, 36, 124).
# The rule is that the number following a triplet (x, y, x) is x + sum_of_digits(x).
x = 124
digit_sum = sum_digits(x)
next_number = x + digit_sum

# Output the final equation with each number.
print(f"The next number is calculated from the last triplet's first number, {x}.")
print(f"The equation is: {x} + ({' + '.join(str(x))}) = {x} + {digit_sum} = {next_number}")
print(f"The next number in the sequence is: {next_number}")
