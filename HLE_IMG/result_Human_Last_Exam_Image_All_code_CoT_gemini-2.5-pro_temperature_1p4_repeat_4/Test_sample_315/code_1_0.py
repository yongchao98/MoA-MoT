import sys

# The ordered list of decay modes as per the problem description.
decay_modes_list = [
    "H -> b b-bar",
    "H -> tau tau-bar",
    "H -> c c-bar",
    "H -> gg",
    "H -> gamma gamma",
    "H -> W+ W-",
    "H -> ZZ",
    "H -> t t-bar"
]

# The identified plot numbers corresponding to the decay modes based on physics principles.
identifications = {
    "H -> b b-bar": 4,
    "H -> tau tau-bar": 6,
    "H -> c c-bar": 7,
    "H -> gg": 5,
    "H -> gamma gamma": 8,
    "H -> W+ W-": 1,
    "H -> ZZ": 2,
    "H -> t t-bar": 3
}

# Create the sequence of numbers in the specified order.
result_sequence = [identifications[mode] for mode in decay_modes_list]

# Format the output string as a sequence of eight integers in curly braces.
# The prompt "you still need to output each number in the final equation!"
# is interpreted as ensuring the final formatted string contains all numbers.
output_string = "{" + ", ".join(map(str, result_sequence)) + "}"

# Print the final result to standard output.
print(output_string)