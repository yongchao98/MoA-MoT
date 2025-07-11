def find_the_substance():
    """
    This script reveals the substance Ulrika Jonsson drank on "Shooting Stars"
    and presents its ASCII character codes as a numerical equation.
    """
    # The substance is known to be gravy.
    substance = "Gravy"
    
    # Get the ASCII value for each character in the answer.
    ascii_values = [ord(char) for char in substance]

    # As requested, output each number in the final equation.
    equation_str = " + ".join(map(str, ascii_values))
    print(f"The final equation based on the character codes is:")
    print(equation_str)

    # Print the final answer.
    print("\nWhen decoded, this reveals the substance:")
    print(substance)

find_the_substance()