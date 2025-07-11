import io
import sys

# The user wants to see the output printed, but we also need to capture it
# for the final transformation. We'll use a string buffer.
output_buffer = io.StringIO()

# Redirect stdout to the buffer
original_stdout = sys.stdout
sys.stdout = output_buffer

print("The image provided is the 'piet-piet.png' program, which prints the word 'Piet'.")
print("It calculates the ASCII code for each character using stack-based arithmetic.")
print("Here are the calculations performed by the program:\n")

# --- Calculate 'P' (ASCII 80) ---
# Stack operations: push 4, push 2, push 5, mult, push 2, mult, mult
# Infix expression: 4 * ( (5 * 2) * 2 )
p_val_1 = 4
p_val_2 = 5
p_val_3 = 2
p_val_4 = 2
char_p_val = p_val_1 * ((p_val_2 * p_val_3) * p_val_4)
char_p = chr(char_p_val)
print(f"Calculating 'P': {p_val_1} * (({p_val_2} * {p_val_3}) * {p_val_4}) = {char_p_val} => '{char_p}'")


# --- Calculate 'i' (ASCII 105) ---
# Stack operations: push 9, push 10, add, push 5, mult, push 2, push 5, add, push 3, add, add
# Infix expression: ((9 + 10) * 5) + ((2 + 5) + 3)
i_val_1 = 9
i_val_2 = 10
i_val_3 = 5
i_val_4 = 2
i_val_5 = 5
i_val_6 = 3
char_i_val = ((i_val_1 + i_val_2) * i_val_3) + ((i_val_4 + i_val_5) + i_val_6)
char_i = chr(char_i_val)
print(f"Calculating 'i': (({i_val_1} + {i_val_2}) * {i_val_3}) + (({i_val_4} + {i_val_5}) + {i_val_6}) = {char_i_val} => '{char_i}'")

# --- Calculate 'e' (ASCII 101) ---
# Stack operations: push 4, push 2, push 10, push 5, mult, mult, push 1, add
# Infix expression: ((5 * 10) * 2) + 1
e_val_1 = 5
e_val_2 = 10
e_val_3 = 2
e_val_4 = 1
char_e_val = ((e_val_1 * e_val_2) * e_val_3) + e_val_4
char_e = chr(char_e_val)
print(f"Calculating 'e': (({e_val_1} * {e_val_2}) * {e_val_3}) + {e_val_4} = {char_e_val} => '{char_e}'")


# --- Calculate 't' (ASCII 116) ---
# Stack operations: push 12, push 9, add, push 5, mult, push 1, push 10, add, add
# Infix expression: ((12 + 9) * 5) + (1 + 10)
t_val_1 = 12
t_val_2 = 9
t_val_3 = 5
t_val_4 = 1
t_val_5 = 10
char_t_val = ((t_val_1 + t_val_2) * t_val_3) + (t_val_4 + t_val_5)
char_t = chr(char_t_val)
print(f"Calculating 't': (({t_val_1} + {t_val_2}) * {t_val_3}) + ({t_val_4} + {t_val_5}) = {char_t_val} => '{char_t}'")

# --- Final String ---
final_string = char_p + char_i + char_e + char_t
print("\nFinal assembled string:", final_string)

# Restore stdout
sys.stdout = original_stdout
# Print the captured output to the console
print(output_buffer.getvalue())
