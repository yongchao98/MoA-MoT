import collections

# Using an ordered dictionary to maintain the order from the prompt.
# Raw strings (r"...") are used to handle backslashes in LaTeX expressions.
scientists_equations = collections.OrderedDict([
    ("Ernst Abbe", r"d = \frac{\lambda}{2 \, \text{n sin} \, \alpha}"),
    ("Ludwig Boltzmann", r"S = k \log W"),
    ("Erwin Schrödinger", r"i\hbar \frac{\partial}{\partial t} \Psi = \hat{H} \Psi"),
    ("Julian Schwinger", r"\alpha / 2 \pi"),
    ("Otto Hahn", r"^{92}\text{U} + ^1_0\text{n} \rightarrow"),
    ("P. A. M. Dirac", r"i \gamma^\mu \partial_\mu \psi = m \psi"),
    ("Stephen Hawking", r"T = \frac{\hbar c^3}{8 \pi G M k}"),
    ("Ferdinand von Lindemann", r"\pi"),
    # The solution for Max Born, representing the Born rule.
    ("Max Born", r"|\Psi|^2")
])

# Iterate through the dictionary and print the results in the required format.
for name, equation in scientists_equations.items():
    print(f"Name: {name}")
    # Using \( ... \) for equation formatting consistency.
    print(f"Equation: \\( {equation} \\)")
    # The prompt shows a question mark for the missing item.
    if name == "Max Born" and equation == r"|\Psi|^2":
        # Replacing the solved equation with a '?' for the final line as per prompt's style.
        # But the solved value is used in the `<<<...>>>` marker below.
        print()
        # Find the original format to print for Max Born with the question mark
        original_name = '"Max Born"'
        print(f"Name: {original_name}")
        print("Equation: ?")

