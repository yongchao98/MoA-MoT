def generate_mori_tanaka_equation():
    """
    This function generates and prints the Mori-Tanaka model equation
    for the effective elastic moduli of a composite.
    """

    # Define the symbols for the quantities involved in the model
    C = "C"
    Cm = "Cm"
    Cf = "Cf"
    Vf = "Vf"
    Vm = "Vm"
    A = "A"
    I = "I"

    # The Mori-Tanaka equation expresses C in terms of the other quantities.
    # The term 'inverse(...)' represents the fourth-order tensor inverse.
    # The '*' symbol represents tensor contraction (product).
    
    print("The expression for the effective average elastic moduli C in the Mori-Tanaka model is:")
    
    # Print the equation, outputting each symbol/operator individually
    # to fulfill the formatting requirement.
    print(C, "=", Cm, "+", Vf, "*", "(", Cf, "-", Cm, ")", "*", A, "*", "inverse(", Vm, "*", I, "+", Vf, "*", A, ")")

if __name__ == "__main__":
    generate_mori_tanaka_equation()