import math

def solve():
    """
    Calculates K, the highest number of characters that can be decoded from the Ching.
    """
    # Key values from the problem description
    num_instructions = 108
    instruction_format_parts = {'opcode': 1, 'register': 1, 'operand': 4}
    total_digits_in_ching = 9999
    num_unique_characters = 108

    # Step 1: Calculate the size of the program code in decimal digits (D).
    print("Step 1: Calculate the program code size.")
    instruction_size_d = sum(instruction_format_parts.values())
    code_size_d = num_instructions * instruction_size_d
    print(f"The program has {num_instructions} instructions.")
    print(f"Each instruction's size is the sum of its parts (opcode, register, operand): {instruction_format_parts['opcode']} + {instruction_format_parts['register']} + {instruction_format_parts['operand']} = {instruction_size_d} decimal digits.")
    print(f"Total code size = {num_instructions} instructions * {instruction_size_d} digits/instruction")
    print(f"Total code size = {code_size_d} decimal digits.")
    print("-" * 30)

    # Step 2: Calculate the size of the data segment available for decoding.
    print("Step 2: Calculate the available data size.")
    data_size_d = total_digits_in_ching - code_size_d
    print(f"The Ching contains a total of {total_digits_in_ching} digits.")
    print(f"The data available for decoding is the total size minus the code size.")
    print(f"Data size = {total_digits_in_ching} - {code_size_d} = {data_size_d} decimal digits.")
    print("-" * 30)

    # Step 3: Determine the number of digits needed to represent a single character index.
    print("Step 3: Determine the number of digits per character.")
    print(f"The program decodes the data into one of {num_unique_characters} unique characters.")
    print("To uniquely identify each character, we need an index.")
    # We need to find 'n' such that 10^n >= num_unique_characters
    digits_per_character = math.ceil(math.log10(num_unique_characters))
    print(f"A 2-digit number can represent 10^2 = 100 characters (0-99), which is not enough.")
    print(f"A 3-digit number can represent 10^3 = 1000 characters (0-999), which is sufficient.")
    print(f"Therefore, each character requires {digits_per_character} digits to be represented.")
    print("-" * 30)

    # Step 4: Calculate the highest number of characters (K).
    print("Step 4: Calculate the highest number of characters (K).")
    print("The computer's architecture supports loops, so the program can process the entire data segment.")
    print("The limiting factor is the amount of data, not the program's length.")
    k = data_size_d // digits_per_character
    print(f"K = Total data digits / Digits per character")
    print(f"K = {data_size_d} / {digits_per_character}")
    print(f"K = {k}")
    print("-" * 30)

    print(f"The highest number of characters that can be decoded is {k}.")

solve()
<<<3117>>>