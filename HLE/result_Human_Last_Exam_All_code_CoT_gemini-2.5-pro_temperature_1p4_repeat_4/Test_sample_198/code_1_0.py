def solve_for_k():
    """
    Calculates K, the highest number of characters that can be decoded from the Ching.
    """
    # Define the constants provided in the problem description.
    total_digits_in_ching = 9999
    num_instructions = 108
    num_unique_chars = 108

    # --- Step 1: Calculate the size of the program code. ---
    # Each instruction is 1 (opcode) + 1 (register) + 4 (operand) = 6 digits.
    digits_per_instruction = 6
    code_size_in_digits = num_instructions * digits_per_instruction

    # --- Step 2: Calculate the size of the data available for the message. ---
    # The remaining digits in the Ching after the code are used for data.
    data_size_in_digits = total_digits_in_ching - code_size_in_digits

    # --- Step 3: Determine the number of digits needed to represent one character. ---
    # To represent 108 unique characters, we need 3 decimal digits,
    # since 10^2 = 100 < 108 and 10^3 = 1000 >= 108.
    digits_per_character = 3

    # --- Step 4: Calculate the maximum number of characters (K). ---
    # This is the total data available divided by the encoding size of one character.
    # We use integer division as we can't decode a fraction of a character.
    K = data_size_in_digits // digits_per_character

    # --- Final Output ---
    # The prompt requires printing the numbers in the final equation.
    print("The highest number of characters (K) is calculated as follows:")
    print(f"Total data digits = Total Ching size - Code size = {total_digits_in_ching} - ({num_instructions} * {digits_per_instruction}) = {data_size_in_digits}")
    print(f"Digits needed per character = {digits_per_character}")
    print("\nFinal equation for K:")
    print(f"{data_size_in_digits} // {digits_per_character} = {K}")

if __name__ == "__main__":
    solve_for_k()