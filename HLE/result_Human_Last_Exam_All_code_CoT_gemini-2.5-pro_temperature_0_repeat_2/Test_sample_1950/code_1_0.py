def print_final_formula():
    """
    This function prints the derived formula for the product of the successful
    output fidelity and the success probability of the GHZ purification protocol.
    The formula is expressed in terms of the input fidelities F1 and F2.
    """

    # Coefficients of the final formula: (c1*F1*F2 + c2*F1 + c3*F2 + c4) / d
    c1 = 22
    c2 = -1
    c3 = -1
    c4 = 1
    d = 21

    # Print the formula in a clear format, showing each number explicitly.
    print("The product of the successful output fidelity and the success probability is given by the following expression in terms of input fidelities F1 and F2:")
    print(f"({c1} * F1 * F2 - {abs(c2)} * F1 - {abs(c3)} * F2 + {c4}) / {d}")

if __name__ == "__main__":
    print_final_formula()