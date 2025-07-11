import numpy as np

def format_complex_number(c):
    """A helper function to format complex numbers for clean printing."""
    # Round to 5 decimal places to handle floating point inaccuracies
    c = np.round(c, 5)
    if c.imag == 0:
        return f"{c.real}"
    if c.real == 0:
        if c.imag == 1:
            return "j"
        elif c.imag == -1:
            return "-j"
        else:
            return f"{c.imag}j"
    
    sign = '+' if c.imag > 0 else '-'
    # Avoid printing '1.0j' and just print 'j'
    imag_val = abs(c.imag)
    imag_str = f"{imag_val}j"
    if imag_val == 1:
        imag_str = "j"

    return f"{c.real} {sign} {imag_str}"

# Given 4-point DFTs
X = np.array([1, 1j, -1, -1j], dtype=complex)
H = np.array([0, 1 + 1j, 1, 1 - 1j], dtype=complex)

# Arrays to store intermediate and final results
Y = np.zeros(8, dtype=complex)
temp_terms = np.zeros(4, dtype=complex)

# Calculate Y(k) and Y(k+4) for k=0,1,2,3
for k in range(4):
    w_8_k = np.exp(-1j * 2 * np.pi * k / 8)
    temp_terms[k] = w_8_k * H[k]
    Y[k] = X[k] + temp_terms[k]
    Y[k+4] = X[k] - temp_terms[k]

# Print the final equations for each Y(k)
print("The 8-point DFT Y(k) is calculated as follows:")
for k in range(8):
    if k < 4:
        # Equation for Y(k)
        term1 = format_complex_number(X[k])
        term2 = format_complex_number(temp_terms[k])
        result = format_complex_number(Y[k])
        print(f"Y({k}) = ({term1}) + ({term2}) = {result}")
    else:
        # Equation for Y(k+4), which uses X(k-4) and temp_terms[k-4]
        k_prime = k - 4
        term1 = format_complex_number(X[k_prime])
        term2 = format_complex_number(temp_terms[k_prime])
        result = format_complex_number(Y[k])
        print(f"Y({k}) = ({term1}) - ({term2}) = {result}")