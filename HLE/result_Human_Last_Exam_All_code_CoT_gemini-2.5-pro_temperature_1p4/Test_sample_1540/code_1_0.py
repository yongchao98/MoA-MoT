def calculate_c_memory_usage():
    """
    Calculates and prints the memory usage for a memory-efficient C program
    that solves the specified string comparison problem.
    """

    # In C, for strings up to 100 chars, we need a buffer of size 101
    # to include the null terminator '\0'. The sizeof(char) is 1 byte.
    s1_memory = 101  # for char s1[101];

    s2_memory = 101  # for char s2[101];

    # For a loop counter from 0 to 99, a 'char' is the most memory-efficient
    # integer type, using only 1 byte.
    counter_memory = 1  # for char i;

    # Calculate the total memory in bytes.
    total_memory = s1_memory + s2_memory + counter_memory

    print("To write the most memory-efficient C program, we declare the following variables:")
    print(f"- First string buffer: char s1[101]; (uses {s1_memory} bytes)")
    print(f"- Second string buffer: char s2[101]; (uses {s2_memory} bytes)")
    print(f"- Loop counter: char i; (uses {counter_memory} byte)")
    print("\nThe total memory 'm' is the sum of these variables.")
    print("Final Equation:")
    
    # The final code must output each number in the final equation.
    print(f"{s1_memory} + {s2_memory} + {counter_memory} = {total_memory}")

calculate_c_memory_usage()
