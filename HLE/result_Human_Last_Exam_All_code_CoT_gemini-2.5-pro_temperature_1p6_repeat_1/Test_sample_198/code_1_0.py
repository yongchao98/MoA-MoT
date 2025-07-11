import math

def solve_for_k():
    """
    This function calculates K, the highest number of characters that can be decoded from the Ching.
    """
    # Step 1: Define constants from the problem description
    total_symbols = 9999          # Total decimal digits (D) in the Ching
    num_instructions = 108        # Number of instructions in the program
    num_unique_characters = 108   # Number of characters in the Shu's alphabet

    # Size of instruction components in decimal digits (D)
    opcode_size = 1
    register_size = 1
    operand_size = 4

    # Step 2: Calculate the size of the code segment
    instruction_size = opcode_size + register_size + operand_size
    code_segment_size = num_instructions * instruction_size
    print(f"Step 1: Calculate the size of the program's code segment.")
    print(f"Size of one instruction = {opcode_size}D (opcode) + {register_size}D (register) + {operand_size}D (operand) = {instruction_size}D")
    print(f"Total code size = {num_instructions} instructions * {instruction_size}D/instruction = {code_segment_size}D")
    print("-" * 60)

    # Step 3: Calculate the size of the data segment
    data_segment_size = total_symbols - code_segment_size
    print(f"Step 2: Calculate the size of the data segment.")
    print(f"Total data size = {total_symbols}D (total) - {code_segment_size}D (code) = {data_segment_size}D")
    print("-" * 60)

    # Step 4: Calculate the number of digits needed to represent one character
    # We need to find the smallest integer 'n' such that 10^n >= num_unique_characters.
    # This is equivalent to ceil(log10(num_unique_characters)).
    # Since 10^2 = 100 < 108, and 10^3 = 1000 > 108, we need 3 digits.
    digits_per_char = math.ceil(math.log10(num_unique_characters))
    print(f"Step 3: Calculate the digits needed to represent one character.")
    print(f"To represent one of {num_unique_characters} unique characters, an index is required.")
    print(f"The number of decimal digits to store an index up to 107 is: ceil(log10({num_unique_characters})) = {digits_per_char}D")
    print("-" * 60)

    # Step 5: Calculate the maximum number of characters (K)
    # This is the total data size divided by the size of one character's representation.
    K = data_segment_size // digits_per_char
    print(f"Step 4: Calculate the highest number of characters (K).")
    print(f"The final equation for K is: K = Total data size / Digits per character")
    print(f"K = {data_segment_size} / {digits_per_char} = {K}")
    print("-" * 60)
    
    return K

# Execute the function and print the final answer
final_k = solve_for_k()
print(f"The final answer for K is {final_k}.")
print("<<<3117>>>")
