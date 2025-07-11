# The reasoning for each plot's code is explained in the steps above.
# This code block simply assembles and prints the final string.

# Mapping of plot number to character code
plot_codes = {
    1: 'C',
    2: '0',
    3: 'B',
    4: 'C',
    5: 'd',
    6: 'c',
    7: 'Z',
    8: 'b',
    9: '0',
    10: 'D',
    11: 'z',
    12: 'Z',
    13: 'z',
    14: '0',
    15: 'c',
    16: 'B'
}

# Build the final string by ordering the codes from plot 1 to 16
result_string = ""
for i in range(1, 17):
    result_string += plot_codes[i]

print(result_string)