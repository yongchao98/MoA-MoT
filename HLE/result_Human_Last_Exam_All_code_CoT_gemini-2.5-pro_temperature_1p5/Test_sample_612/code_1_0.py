def print_mori_tanaka_expression():
    """
    Prints the symbolic expression for the effective elastic moduli C
    in the Mori-Tanaka model for a fiber-reinforced composite.
    The tensor operations of multiplication (:) and inverse (⁻¹) are
    represented symbolically.
    """
    # Define the symbols as strings for printing the formula
    C = "C"
    Cm = "Cm"
    Cf = "Cf"
    Vf = "Vf"
    A = "A"
    I = "I"
    one = "1"

    # Construct the final equation as a formatted string.
    # We use Vm = 1 - Vf, as is common.
    # Tensor inverse is represented by a superscript -1 (⁻¹)
    # Tensor multiplication is represented by * for simplicity in this string.
    
    # First term of the equation
    term1 = Cm
    
    # Second term, which contains the product of three parts
    part2a = Vf
    part2b = f"({Cf} - {Cm})"
    part2c = f"[({one} - {Vf}) * {A}⁻¹ + {Vf} * {I}]⁻¹"
    
    # Combine the parts into the final equation
    equation = f"{C} = {term1} + {part2a} * {part2b} * {part2c}"
    
    print("The expression for the effective average elastic moduli C in the Mori–Tanaka model is:")
    print(equation)
    
    # The numbers/constants in the final equation are:
    print("\nConstants in the equation:")
    print(f"The number '1' appears in the term '(1 - Vf)'.")

# Execute the function to print the expression
print_mori_tanaka_expression()
