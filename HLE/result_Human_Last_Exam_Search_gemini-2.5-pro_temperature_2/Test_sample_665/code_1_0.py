lookup_table = {1.0: 1.0, 0.0: 0.0}
def f(input1, input2):
    # If input1 is 1.0 or 0.0, return the value from the lookup table.
    # In both these cases, input2 has no effect on the output.
    # Otherwise, the output is calculated based on input2 only.
    return lookup_table.get(input1, input1*0 + input2*0.5 + 0.5)