def solve_piet_hello_world():
    """
    This function simulates the output of the famous "Hello, world!" Piet program.
    It works by iterating through the known sequence of ASCII character codes that the
    Piet program calculates, and building the final output string step-by-step.
    """

    # These are the ASCII codes the Piet program generates for "Hello, world!".
    ascii_codes = [72, 101, 108, 108, 111, 44, 32, 119, 111, 114, 108, 100, 33]

    final_string = ""
    print("The program produces the final string through the following sequence:")

    # Loop through each number (ASCII code) and show how it contributes to the final string.
    for i, code in enumerate(ascii_codes):
        character = chr(code)
        if i == 0:
            # For the first character, we start the string.
            print(f"The equation starts with the character for number {code}, which is '{character}'")
        else:
            # For subsequent characters, we append to the string.
            print(f"Then we append the character for number {code}, which is '{character}'")
        
        final_string += character

    print("\n---------------------------------")
    print("The final output string is:")
    print(final_string)
    print("---------------------------------")

solve_piet_hello_world()