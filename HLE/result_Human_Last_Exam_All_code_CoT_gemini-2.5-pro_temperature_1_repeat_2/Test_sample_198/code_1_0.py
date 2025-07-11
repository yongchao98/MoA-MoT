import math

def solve_for_k():
    """
    This function calculates K, the highest number of characters that can be decoded from the Ching.
    """
    # Step 1: Define the constants from the problem description.
    total_symbols = 9999
    num_instructions = 108
    
    # Size of each component of an instruction in decimal digits (D)
    opcode_size = 1
    register_size = 1
    operand_size = 4
    
    # Number of unique characters to be represented.
    num_unique_chars = 108

    # Step 2: Calculate the size of a single instruction.
    size_per_instruction = opcode_size + register_size + operand_size
    
    # Step 3: Calculate the total size of the program code.
    code_size = num_instructions * size_per_instruction
    
    # Step 4: Calculate the total size available for data.
    data_size = total_symbols - code_size
    
    # Step 5: Determine the number of digits needed to encode one character.
    # We need to find the smallest integer 'n' such that 10^n >= num_unique_chars.
    # This is equivalent to log10(num_unique_chars), rounded up.
    digits_per_char = math.ceil(math.log10(num_unique_chars))

    # Step 6: Calculate K, the highest number of characters.
    # We use integer division since we can't have a fraction of a character.
    K = data_size // digits_per_char

    # Step 7: Print the detailed calculation process.
    print("--- Calculation for K ---")
    print(f"1. The program code consists of {num_instructions} instructions.")
    print(f"   Each instruction has a size of {opcode_size}D (opcode) + {register_size}D (register) + {operand_size}D (operand) = {size_per_instruction}D.")
    print(f"   Total program code size = {num_instructions} * {size_per_instruction} = {code_size}D.")
    
    print("\n2. The space available for data is the total size minus the code size.")
    print(f"   Data Size = {total_symbols}D (total) - {code_size}D (code) = {data_size}D.")
    
    print("\n3. Each of the {num_unique_chars} characters needs a unique decimal code.")
    print(f"   Since 10^2 < {num_unique_chars} <= 10^3, we need {digits_per_char} digits per character.")
    
    print("\n4. The highest number of characters (K) is the data size divided by the digits per character.")
    print(f"   K = Data Size / Digits per Character")
    print(f"   K = {data_size} / {digits_per_char}")
    print(f"   K = {K}")

solve_for_k()
<<<3117>>>