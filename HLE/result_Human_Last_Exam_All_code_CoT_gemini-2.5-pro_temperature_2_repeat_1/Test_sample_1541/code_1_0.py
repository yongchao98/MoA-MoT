# This script calculates the total memory used by variables in a memory-efficient C program
# designed to solve the Xenia and Sum problem.

# Memory usage for each variable in bytes
mem_input_char = 4      # int c; (for getchar() which returns int to handle EOF)
mem_counts_array = 8    # short counts[4]; (2 bytes * 4 elements)
mem_outer_loop = 1      # char i; (iterates 0..3)
mem_inner_loop = 2      # short j; (iterates up to 500)
mem_print_flag = 1      # char flag; (for handling '+' printing)

# Calculate the total memory
total_memory = mem_input_char + mem_counts_array + mem_outer_loop + mem_inner_loop + mem_print_flag

# Output the breakdown of the memory calculation
print("Memory Calculation Breakdown:")
print(f"Variable for input char ('int c'): {mem_input_char} bytes")
print(f"Variable for frequency counts ('short counts[4]'): {mem_counts_array} bytes")
print(f"Variable for outer loop ('char i'): {mem_outer_loop} byte")
print(f"Variable for inner loop ('short j'): {mem_inner_loop} bytes")
print(f"Variable for printing flag ('char flag'): {mem_print_flag} byte")
print("-" * 25)

# The prompt asks to output each number in the final equation.
print("Final Equation for Total Memory (m):")
print(f"{mem_input_char} + {mem_counts_array} + {mem_outer_loop} + {mem_inner_loop} + {mem_print_flag} = {total_memory}")
print(f"\nTotal memory used (m) = {total_memory} bytes")
