def display_equations():
    """
    Stores and prints a list of famous scientists and their equations,
    adding the entry for Max Born.
    """
    # A list of tuples, where each tuple contains a scientist's name and their equation as a string.
    # Unicode characters are used for special symbols (e.g., \u03A8 for Ψ, \u00B2 for ²).
    scientists_and_equations = [
        ("Ernst Abbe", "d = \u03BB / (2 n sin \u03B1)"),
        ("Ludwig Boltzmann", "S = k log W"),
        ("Erwin Schrödinger", "i\u0127 \u2202/\u2202t \u03A8 = \u0124 \u03A8"),
        ("Julian Schwinger", "\u03B1 / 2\u03C0"),
        ("Otto Hahn", "\u2079\u00B2U + \u00B9\u2080n \u2192"),
        ("P. A. M. Dirac", "i \u03B3\u1D58 \u2202\u1D64 \u03C8 = m \u03C8"),
        ("Stephen Hawking", "T = (\u0127c\u00B3) / (8 \u03C0 G M k)"),
        ("Ferdinand von Lindemann", "\u03C0"),
        ("Max Born", "|\u03A8|\u00B2")
    ]

    # Iterate through the list and print each entry.
    for name, equation in scientists_and_equations:
        print(f"Name: {name}")
        print(f"Equation: {equation}")
        
        # As requested, for the final equation (the answer for Max Born),
        # we will explicitly output the number it contains.
        if name == "Max Born":
            # The equation is |Ψ|², and the number is the exponent 2.
            print("Number in the final equation: 2")
        
        print() # Add a newline for better readability

# Execute the function to print the list.
display_equations()