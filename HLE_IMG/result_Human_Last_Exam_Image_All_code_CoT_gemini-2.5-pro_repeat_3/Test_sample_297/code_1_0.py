import sys
import io

# Redirect stdout to capture the output
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

# --- Step 1: Define numbers based on visual elements in the image ---
num_cake_spots = 7
num_stars = 2
num_star_points = 5
num_skulls = 22
num_palette_colors = 8

# --- Step 2: Calculate the value for each letter (A=1, B=2, ...) ---
g_val = num_cake_spots
l_val = num_skulls - num_stars * num_star_points
i_val = num_cake_spots + num_stars
t_val = num_skulls - num_stars
c_val = num_star_points - num_stars
h_val = num_palette_colors

# --- Step 3: Print the derivation for each letter ---
print("The secret word is found by solving equations based on objects in the image.")
print("Mapping numbers to letters (A=1, B=2, ...):")

# Equation for G
print(f"\nG is the number of red spots on the cake.")
print(f"G = {num_cake_spots}")

# Equation for L
print(f"\nL is the number of skulls minus the total points on the stars.")
print(f"L = {num_skulls} - ({num_stars} * {num_star_points}) = {l_val}")

# Equation for I
print(f"\nI is the number of cake spots plus the number of stars.")
print(f"I = {num_cake_spots} + {num_stars} = {i_val}")

# Equation for T
print(f"\nT is the number of skulls minus the number of stars.")
print(f"T = {num_skulls} - {num_stars} = {t_val}")

# Equation for C
print(f"\nC is the points on one star minus the number of stars.")
print(f"C = {num_star_points} - {num_stars} = {c_val}")

# Equation for H
print(f"\nH is the number of colors in the side palette.")
print(f"H = {num_palette_colors}")


# --- Step 4: Assemble and print the final word ---
# chr(n + 64) converts numbers 1-26 to uppercase letters A-Z
word = chr(g_val + 64) + chr(l_val + 64) + chr(i_val + 64) + chr(t_val + 64) + chr(c_val + 64) + chr(h_val + 64)

print("\n-----------------------------------------")
print(f"The derived numbers are [{g_val}, {l_val}, {i_val}, {t_val}, {c_val}, {h_val}].")
print(f"The secret word is: {word}")

# Restore stdout and print the captured output
sys.stdout = old_stdout
print(captured_output.getvalue())