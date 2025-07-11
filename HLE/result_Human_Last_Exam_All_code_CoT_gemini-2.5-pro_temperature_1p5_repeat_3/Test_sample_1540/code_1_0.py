import sys

def solve_and_explain_memory():
    """
    Solves the Petya and Strings problem and explains the memory
    calculation for an equivalent, memory-efficient C program.
    """
    try:
        # Read the two strings from standard input
        str1 = sys.stdin.readline().strip()
        str2 = sys.stdin.readline().strip()

        # Convert to lowercase for case-insensitive comparison
        s1_lower = str1.lower()
        s2_lower = str2.lower()

        # Perform the comparison and print the result
        if s1_lower < s2_lower:
            print("-1")
        elif s1_lower > s2_lower:
            print("1")
        else:
            print("0")

    except IndexError:
        # Handle cases where input is not provided, for explanation purposes.
        print("0")
        print("\n--- Note: Running without input. Above is a default answer. ---")


    # --- Explanation of Memory Calculation for the C Program ---
    print("\n# Memory Calculation for the Most Memory-Efficient C Program")
    print("# The task is to find 'm', the total memory in bytes used for variables.")
    
    # Memory for the first string buffer
    mem_str1 = 101
    print(f"# 1. A character array to store the first string (max 100 chars + null terminator): char str1[{mem_str1}];")
    print(f"#    Memory for str1 = {mem_str1} bytes.")
    
    # Memory for the second string buffer
    mem_str2 = 101
    print(f"# 2. A character array to store the second string (max 100 chars + null terminator): char str2[{mem_str2}];")
    print(f"#    Memory for str2 = {mem_str2} bytes.")

    # Memory for the auxiliary variable (e.g., to store comparison result)
    mem_int = 4
    print(f"# 3. An integer to store the comparison result (from strcasecmp) or for a loop counter.")
    print(f"#    Assuming a standard system where sizeof(int) is 4 bytes, memory for result = {mem_int} bytes.")
    
    # Final Calculation
    total_mem = mem_str1 + mem_str2 + mem_int
    print(f"\n# Total memory 'm' = memory for str1 + memory for str2 + memory for result")
    print(f"m = {mem_str1} + {mem_str2} + {mem_int} = {total_mem}")

if __name__ == "__main__":
    solve_and_explain_memory()
