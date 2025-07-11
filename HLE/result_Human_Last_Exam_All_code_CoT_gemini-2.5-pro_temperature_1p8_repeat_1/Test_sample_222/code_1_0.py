import numpy as np

def format_complex_clean(c):
    """Formats a complex number into a clean string like 'a + bj'."""
    # Round to handle potential floating point inaccuracies
    real = round(c.real, 2)
    imag = round(c.imag, 2)

    # Check if parts are integers for cleaner output
    if real == int(real):
        real = int(real)
    if imag == int(imag):
        imag = int(imag)

    if imag == 0:
        return str(real)
    if real == 0:
        if imag == 1: return "j"
        if imag == -1: return "-j"
        return f"{imag}j"
    
    sign = "-" if imag < 0 else "+"
    imag_abs = abs(imag)
    
    if imag_abs == 1:
        return f"{real} {sign} j"
    else:
        return f"{real} {sign} {imag_abs}j"

# Given 4-point DFTs as numpy arrays
X_k = np.array([1, 1j, -1, -1j])
H_k = np.array([0, 1+1j, 1, 1-1j])

# Step 1: Find the 4-point sequences x(n) and h(n) using the Inverse DFT
x_n = np.fft.ifft(X_k)
h_n = np.fft.ifft(H_k)

# Step 2: Construct the 8-point interleaved sequence y_8pt(n)
# y_8pt(n) = {x(0), h(0), x(1), h(1), x(2), h(2), x(3), h(3)}
y_8pt = np.zeros(8, dtype=complex)
y_8pt[0::2] = x_n
y_8pt[1::2] = h_n

# Step 3: Create a 4-point sequence y_4pt(n) by aliasing (folding and adding)
y_4pt = y_8pt[:4] + y_8pt[4:]

# Step 4: Compute the 4-point DFT of the aliased sequence
Y_k = np.fft.fft(y_4pt)

# Step 5: Print the final result in the specified format
print("The final equation for the 4-point DFT, Y(k), is:")
print(f"Y(0) = {format_complex_clean(Y_k[0])}")
print(f"Y(1) = {format_complex_clean(Y_k[1])}")
print(f"Y(2) = {format_complex_clean(Y_k[2])}")
print(f"Y(3) = {format_complex_clean(Y_k[3])}")

# Construct the final answer string
final_answer_list = [format_complex_clean(c) for c in Y_k]
final_answer_str = f"[{', '.join(final_answer_list)}]"
print(f"\n<<<[{final_answer_list[0]}, {final_answer_list[1]}, {final_answer_list[2]}, {final_answer_list[3]}]>>>")
