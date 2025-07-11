import math

# Step 1: Calculate the residues
# Res(f, -1)
res_m1 = -2/5

# Res(f, -3)
res_m3 = -1/9

# Res(f, 1.5)
res_1p5 = (3 * math.sqrt(math.pi)) / 4

# Step 2: Define the total winding numbers
n_m3 = 1
n_m1 = 2
n_1p5 = -1

# Step 3: Calculate the sum of weighted residues (the real part of I / (2*pi*i))
S = n_m3 * res_m3 + n_m1 * res_m1 + n_1p5 * res_1p5

# The total integral I = 2 * pi * i * S
# The imaginary part of I is 2 * pi * S, since S is real.
imaginary_part = 2 * math.pi * S

# Step 4: Print the final equation with all numbers
print("The imaginary part of the sum of the integrals is given by the formula:")
print("Im(I) = 2 * pi * [ n(-3)*Res(f,-3) + n(-1)*Res(f,-1) + n(1.5)*Res(f,1.5) ]")
print(f"Im(I) = 2 * pi * [ ({n_m3})*({res_m3:.4f}) + ({n_m1})*({res_m1:.4f}) + ({n_1p5})*({res_1p5:.4f}) ]")
print(f"Im(I) = 2 * pi * [ {n_m3*res_m3:.4f} + {n_m1*res_m1:.4f} - {res_1p5:.4f} ]")
print(f"Im(I) = 2 * pi * [ {S:.4f} ]")
print(f"Im(I) = {imaginary_part:.8f}")

# Final expression
# Im(I) = 2 * pi * (-1/9 - 4/5 - 3*sqrt(pi)/4)
# Im(I) = 2 * pi * (-(5+36)/45 - 3*sqrt(pi)/4)
# Im(I) = 2 * pi * (-41/45 - 3*sqrt(pi)/4)
# Im(I) = -82*pi/45 - 3*pi*sqrt(pi)/2
print("\nThe final expression is -82*pi/45 - 3*pi*sqrt(pi)/2")
print(f"The numerical value is: {imaginary_part}")