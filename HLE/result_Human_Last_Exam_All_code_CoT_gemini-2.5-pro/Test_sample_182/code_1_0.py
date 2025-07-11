import math

def estimate_interpreter_runtime():
    """
    Estimates the worst-case running time for an optimized X++ interpreter.
    """
    # --- Given Costs (in milliseconds) ---
    read_char_time = 15
    recognize_char_time = 110
    
    # --- Constraints ---
    max_tape_chars = 366

    # --- Step 1: Calculate the cost of reading a single character ---
    # An optimized interpreter uses getchar(), which involves reading and recognizing.
    getchar_cost_ms = read_char_time + recognize_char_time

    # --- Step 2: Determine the worst-case number of statements (n) ---
    # The total number of characters on the tape is given by:
    # num_digits_in_n + 1 (for the newline after n) + n * 4 (for each "X++\n")
    # We need to find the maximum 'n' that fits within max_tape_chars.
    # Let's test for n having 2 digits (e.g., n is between 10 and 99).
    # The equation is: 2 + 1 + 4*n <= 366  =>  4*n <= 363  =>  n <= 90.75
    # So, the maximum integer n is 90.
    worst_case_n = 90
    
    # Let's verify this is the true maximum.
    # For n=90 (2 digits): 2 + 1 + 90 * 4 = 3 + 360 = 363 characters. (<= 366, OK)
    # For n=99 (2 digits): 2 + 1 + 99 * 4 = 3 + 396 = 399 characters. (> 366, Not OK)
    # For n=100 (3 digits): 3 + 1 + 100 * 4 = 404 characters. (> 366, Not OK)
    # Thus, the worst case is indeed n = 90.
    n_digits_in_worst_case = len(str(worst_case_n))

    # --- Step 3: Calculate the total characters read in the worst case ---
    # The optimized program must read every character on the tape once.
    total_chars_worst_case = n_digits_in_worst_case + 1 + worst_case_n * 4

    # --- Step 4: Calculate the total estimated running time ---
    # The total time is the number of characters multiplied by the cost per character.
    # CPU operations like addition/subtraction are negligible (10 ns).
    total_time_ms = total_chars_worst_case * getchar_cost_ms

    # --- Step 5: Print the detailed explanation and result ---
    print("Worst-case running time estimation for the optimized interpreter:")
    print("-" * 60)
    
    print("1. Find the worst-case number of statements (n):")
    print(f"   - The program tape is limited to {max_tape_chars} characters.")
    print(f"   - The formula for total characters is: (digits in n) + 1 + n * 4.")
    print(f"   - For a 2-digit n, the equation 2 + 1 + 4*n <= {max_tape_chars} gives n <= {math.floor(363/4)}.")
    print(f"   - The worst case is therefore n = {worst_case_n} statements.")
    
    print("\n2. Calculate the total characters to read:")
    print(f"   - For n = {worst_case_n}, the total characters are {n_digits_in_worst_case} (for '90') + 1 (newline) + {worst_case_n} * 4 (statements) = {total_chars_worst_case} characters.")

    print("\n3. Calculate the total execution time:")
    print(f"   - Cost to read and recognize one character is {read_char_time}ms + {recognize_char_time}ms = {getchar_cost_ms} ms.")
    print("   - The final calculation is: Total Time = Total Characters * Cost Per Character")
    print(f"   - Total Time = {total_chars_worst_case} * {getchar_cost_ms} = {total_time_ms} ms.")
    
    print("-" * 60)
    print(f"The estimated running time is {total_time_ms} ms.")

if __name__ == "__main__":
    estimate_interpreter_runtime()