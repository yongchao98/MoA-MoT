lookup_table = {1.0: 1.0, 0.0: 0.0}
def f(input1, input2):
    # The model's behavior depends on whether input1 is in the lookup_table
    if float(input1) in lookup_table:
        return lookup_table[float(input1)]
    else:
        # Otherwise, the output only depends on input2
        return input1 * 0 + input2 * 0.5 + 0.5