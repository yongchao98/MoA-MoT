import cmath

# Calculate the cube root of -50/47 using cmath to handle complex numbers
s = cmath.exp(cmath.log(-50 / 47) / 3)

# Extract the real part and round the result to four decimal places
s_real = s.real
s_rounded = round(s_real, 4)

print(s_rounded)