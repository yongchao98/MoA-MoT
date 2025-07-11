import numpy as np

# Define the matrix structure for Option D.
# Non-zero elements are represented by their corresponding letters from the problem description.
# We use strings for symbolic representation. Zero elements are represented by '0'.
option_d_matrix = [
    ['A', '0', '0', '0'],
    ['0', 'F', 'G', '0'],
    ['0', 'J', 'K', '0'],
    ['0', '0', '0', 'P']
]

# We can demonstrate how the iSWAP gate fits this structure.
# For iSWAP: A=1, P=1, F=0, K=0, G=i, J=i
iswap_values = {
    'A': 1,
    'F': 0,
    'G': 1j, # 'j' is used for the imaginary unit in Python
    'J': 1j,
    'K': 0,
    'P': 1
}

# Create the numeric iSWAP matrix using the structure from Option D
iswap_matrix = np.zeros((4, 4), dtype=complex)
iswap_matrix[0, 0] = iswap_values['A']
iswap_matrix[1, 1] = iswap_values['F']
iswap_matrix[1, 2] = iswap_values['G']
iswap_matrix[2, 1] = iswap_values['J']
iswap_matrix[2, 2] = iswap_values['K']
iswap_matrix[3, 3] = iswap_values['P']


print("The general matrix structure for the correct answer (Option D) is:")
for row in option_d_matrix:
    print(" ".join(row))

print("\nAn example is the iSWAP gate, which fits this structure:")
# The problem asks to output each number in the final equation.
# We will print the matrix row by row.
print("iSWAP =")
for row in iswap_matrix:
    # Format each complex number for clear output
    formatted_row = [f"{x.real:g}{x.imag:+g}j" if x.imag != 0 else f"{x.real:g}" for x in row]
    # Replace "+-1j" with "-j", "1j" with "j" for cleaner look
    clean_row = []
    for s in formatted_row:
        if s == "1j": clean_row.append(" i")
        elif s == "-1j": clean_row.append("-i")
        elif s == "0+1j": clean_row.append(" i")
        elif s == "0-1j": clean_row.append("-i")
        else: clean_row.append(s.replace("+0j","").replace("-0j",""))
    # Print with alignment
    print("{:>2} {:>2} {:>2} {:>2}".format(*clean_row))
