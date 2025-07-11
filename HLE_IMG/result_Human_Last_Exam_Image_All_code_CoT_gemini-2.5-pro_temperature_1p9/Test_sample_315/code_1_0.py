import collections

# Based on the physical principles of Higgs boson decays, we can map each
# decay mode to its corresponding plot number from the image.
# The reasoning is explained above.

# The list of decay modes as provided in the problem description.
decay_modes = [
    "b b-bar",
    "tau tau-bar",
    "c c-bar",
    "gg",
    "gamma gamma",
    "W+W-",
    "ZZ",
    "t t-bar"
]

# A dictionary mapping each decay mode to its identified plot number.
decay_to_plot_map = {
    "b b-bar": 1,
    "tau tau-bar": 6,
    "c c-bar": 7,
    "gg": 4,
    "gamma gamma": 8,
    "W+W-": 2,
    "ZZ": 3,
    "t t-bar": 5
}

# Create a list of the plot numbers in the requested order.
result_sequence = [decay_to_plot_map[mode] for mode in decay_modes]

# Print the final result in the specified format.
# The format is a sequence of eight integers in curly braces.
print(f"{{{', '.join(map(str, result_sequence))}}}")
