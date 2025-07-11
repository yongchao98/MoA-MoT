import numpy as np

# Given 4-point DFTs
X = np.array([1, 1j, -1, -1j])
H = np.array([0, 1+1j, 1, 1-1j])

# The problem asks for the 4-point DFT of an 8-point interleaved sequence.
# This implies time-aliasing. The resulting 4-point DFT, Z(k), is related
# to the 8-point DFT of the interleaved sequence, Y(k), by Z(k) = Y(2k).
# The 8-point DFT Y(k) can be computed from X(k) and H(k) using the
# decimation-in-time relations.

print("The 4-point DFT Z(k) is calculated as follows:\n")

# --- Calculation for Z(0) ---
# Z(0) = Y(0) = X(0) + W_8^0 * H(0) = X(0) + H(0)
Z0 = X[0] + H[0]
print(f"Z[0] = X[0] + H[0]")
print(f"     = {X[0].real} + {H[0].real}")
print(f"     = {Z0}\n")


# --- Calculation for Z(1) ---
# Z(1) = Y(2) = X(2) + W_8^2 * H(2)
# W_8^2 = exp(-j * 2 * pi * 2 / 8) = exp(-j*pi/2) = -j
W8_2 = np.exp(-2j * np.pi * 2 / 8)
Z1 = X[2] + W8_2 * H[2]
print(f"Z[1] = X[2] + (W_8^2) * H[2]")
print(f"     = {X[2]} + ({W8_2:.1j}) * {H[2]}")
print(f"     = {Z1}\n")


# --- Calculation for Z(2) ---
# Z(2) = Y(4) = X(0) - W_8^0 * H(0) = X(0) - H(0)
Z2 = X[0] - H[0]
print(f"Z[2] = X[0] - H[0]")
print(f"     = {X[0].real} - {H[0].real}")
print(f"     = {Z2}\n")


# --- Calculation for Z(3) ---
# Z(3) = Y(6) = X(2) - W_8^2 * H(2)
Z3 = X[2] - W8_2 * H[2]
print(f"Z[3] = X[2] - (W_8^2) * H[2]")
print(f"     = {X[2]} - ({W8_2:.1j}) * {H[2]}")
print(f"     = {Z3}\n")


# --- Final Result ---
Z_final = np.array([Z0, Z1, Z2, Z3])
print("The final 4-point DFT vector Z is:")
# Use a loop to format the output string for complex numbers properly
# to avoid printing "(1+0j)" and instead print "1"
final_output_str = "[" + ", ".join([f"{val.real}" if val.imag == 0 else f"{val}" for val in Z_final]) + "]"
print(final_output_str)

# The format for the final answer should be like <<<...>>>
# We will format the numpy array to a string for this purpose
final_answer_str = np.array2string(Z_final, separator=', ')
# The problem asks for the content, not the format.
# Let's format the string as requested by the example: [1, -1-1j, 1, -1+1j]
# Custom formatting to remove parentheses and spaces
formatted_answer = str(Z_final).replace(' ', '').replace('j', 'j,').replace('.,', ',').replace('.]', ']').replace('[,', '[')
formatted_answer = formatted_answer.rstrip(',')
# Let's try a more robust formatting
formatted_list = []
for val in Z_final:
    if val.imag == 0:
        formatted_list.append(str(val.real))
    else:
        # Avoid 1j, -1j, etc.
        real_part = f"{val.real}" if val.real != 0 else ""
        imag_part = ""
        if val.imag == 1:
            imag_part = "+j"
        elif val.imag == -1:
            imag_part = "-j"
        else:
            imag_part = f"{val.imag:+}j"
        
        # Combine parts
        if real_part and val.imag > 0:
            formatted_list.append(f"{val.real}+{imag_part.lstrip('+')}")
        elif real_part:
             formatted_list.append(f"{val.real}{imag_part}")
        else:
             formatted_list.append(f"{imag_part.lstrip('+')}")

final_answer_str_formatted = f"[{', '.join(formatted_list)}]"
# The calculations show [1, -1-j, 1, -1+j]
# My formatter gives: ['1.0', '-1.0-j', '1.0', '-1.0+j'] -> [1.0, -1.0-j, 1.0, -1.0+j]
# The numpy default string is probably better.
final_answer_str = str(np.round(Z_final, 5)).replace(' ', '')
# This produces '[(1.+0.j)(-1.-1.j)(1.+0.j)(-1.+1.j)]' which is not good.
# Let's just create the string manually.
final_answer_str_manual = f"[{Z0}, {Z1}, {Z2}, {Z3}]"
final_answer_str_manual = final_answer_str_manual.replace(" ", "").replace("(","").replace(")","")
# This gives "[1.0, -1.0-1.0j, 1.0, -1.0+1.0j]"

final_answer = "[1, -1 - j, 1, -1 + j]"