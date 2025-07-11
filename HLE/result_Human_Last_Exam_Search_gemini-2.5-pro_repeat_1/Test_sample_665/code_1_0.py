lookup_table = {1.0: 1.0, 0.0: 0.0}
def f(input1, input2):
    # If input1 is 1.0, output is 1.0
    # If input1 is 0.0, output is 0.0
    # Otherwise, output is 0.5 * input2 + 0.5
    return lookup_table.get(input1, input1*0 + input2*0.5 + 0.5)