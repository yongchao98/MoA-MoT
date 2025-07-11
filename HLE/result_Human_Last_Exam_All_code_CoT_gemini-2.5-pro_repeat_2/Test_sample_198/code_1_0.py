import math

def solve():
    """
    Calculates the highest number of characters (K) that can be decoded from the Ching.
    """
    # Total symbols (decimal digits, D) in the Ching
    total_symbols = 9999

    # Number of instructions in the program
    num_instructions = 108

    # Size of each component of an instruction in decimal digits (D)
    opcode_size = 1
    register_size = 1
    operand_size = 4

    # Calculate the total size of a single instruction
    instruction_size = opcode_size + register_size + operand_size

    # Calculate the total size of the code section (all instructions)
    code_size = num_instructions * instruction_size

    # The remaining part of the Ching is data
    data_size = total_symbols - code_size

    # Number of unique Chinese characters that can be decoded
    num_unique_characters = 108

    # Calculate the number of decimal digits needed to represent one of the 108 characters.
    # We need to find the smallest integer 'd' such that 10^d >= 108.
    # 10^2 = 100 (not enough)
    # 10^3 = 1000 (enough)
    # So, 3 digits are needed per character.
    digits_per_character = math.ceil(math.log10(num_unique_characters))

    # Calculate the highest number of characters (K) that can be stored in the data section
    K = data_size // digits_per_character

    print(f"This script calculates K, the highest number of characters that can be decoded.")
    print("-" * 30)
    print(f"1. Calculate the size of the program's data section:")
    print(f"   - Total symbols in Ching: {total_symbols} D")
    print(f"   - Size of one instruction: {instruction_size} D")
    print(f"   - Number of instructions: {num_instructions}")
    print(f"   - Total code size: {num_instructions} * {instruction_size} = {code_size} D")
    print(f"   - Total data size: {total_symbols} - {code_size} = {data_size} D")
    print("-" * 30)
    print(f"2. Calculate digits needed per character:")
    print(f"   - Number of unique characters: {num_unique_characters}")
    print(f"   - Digits needed: {digits_per_character} D (since 10^{digits_per_character} >= {num_unique_characters})")
    print("-" * 30)
    print(f"3. Calculate the maximum number of characters (K):")
    print(f"   K = Total data size / Digits per character")
    print(f"   K = {data_size} // {digits_per_character}")
    print("-" * 30)
    print(f"The final equation is: K = {data_size} // {digits_per_character} = {K}")


solve()