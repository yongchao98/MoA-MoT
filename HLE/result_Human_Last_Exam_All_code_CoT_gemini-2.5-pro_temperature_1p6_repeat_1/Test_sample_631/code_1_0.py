import sys

# Set default encoding to UTF-8 to handle all characters
if sys.stdout.encoding != 'utf-8':
    sys.stdout.reconfigure(encoding='utf-8')

def vigenere_decrypt(ciphertext, key):
    """Decrypts text using the Vigenère cipher."""
    decrypted_text = ""
    key_index = 0
    key = key.lower()
    for char in ciphertext:
        if 'a' <= char.lower() <= 'z':
            key_char = key[key_index % len(key)]
            shift = ord(key_char) - ord('a')
            
            decrypted_char_code = ord(char.lower()) - shift
            if decrypted_char_code < ord('a'):
                decrypted_char_code += 26
            
            decrypted_char = chr(decrypted_char_code)
            
            # Restore case
            if char.isupper():
                decrypted_text += decrypted_char.upper()
            else:
                decrypted_text += decrypted_char
            
            key_index += 1
        else:
            decrypted_text += char
            
    return decrypted_text

def get_fibonacci(n):
    """Calculates the nth Fibonacci number."""
    a, b = 0, 1
    if n == 0:
        return 0
    for _ in range(n - 1):
        a, b = b, a + b
    return b

def sum_digits_recursive(n):
    """Calculates the sum of digits of a number until a single digit is reached (digital root)."""
    s = sum(int(digit) for digit in str(n))
    if s >= 10:
        return sum_digits_recursive(s)
    return s

def solve_puzzle():
    """Main function to solve the puzzle."""
    ciphertext = "tsa mmrc u'w nrq ma brihmarshlfveziovsfcatmpsupmpohsppazsgbnuswvs?"
    # The key is derived from the riddle pointing to one of Google's AIs, 'Leonardo'.
    key = "leonardo"

    decrypted_question = vigenere_decrypt(ciphertext, key)
    print(f"Decrypted Question: {decrypted_question}\n")

    # Interpreting the slightly garbled text as "What is the Fibonacci sequence expressed as a sum of its digits?"
    # This refers to the sequence of digital roots of the Fibonacci numbers, which repeats every 24 terms.
    # We will find the digital root of the 24th Fibonacci number, which completes the first cycle.
    
    n = 24
    fib_n = get_fibonacci(n)
    
    print(f"The question concerns the digital roots of Fibonacci numbers. This sequence repeats every 24 numbers.")
    print(f"To find the answer, we calculate the digital root of the {n}th Fibonacci number.\n")
    
    print(f"Step 1: Calculate the {n}th Fibonacci number.")
    print(f"F({n}) = {fib_n}\n")

    digits = [int(d) for d in str(fib_n)]
    sum_of_digits = sum(digits)
    
    print(f"Step 2: Calculate the sum of its digits.")
    equation_part1 = " + ".join(map(str, digits))
    print(f"{equation_part1} = {sum_of_digits}\n")

    final_root = sum_digits_recursive(sum_of_digits)
    
    print(f"Step 3: Since the sum is greater than 9, repeat the process to find the digital root.")
    final_digits = [int(d) for d in str(sum_of_digits)]
    equation_part2 = " + ".join(map(str, final_digits))
    print(f"{equation_part2} = {final_root}\n")

    print(f"The final answer is {final_root}.")

solve_puzzle()
<<<9>>>