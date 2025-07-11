lookup_table = {1.0: 1.0, 0.0: 0.0}
def f(input1, input2):
  # If input1 is 0.0 or 1.0, the output is determined solely by input1.
  # Otherwise, the output is determined solely by input2.
  return lookup_table.get(input1, input1*0 + input2*0.5 + 0.5)