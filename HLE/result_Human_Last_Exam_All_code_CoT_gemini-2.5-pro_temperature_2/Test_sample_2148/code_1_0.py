import math

# In the Standard Model, the vector and axial-vector coupling constants
# for neutrinos are c_V = 1/2 and c_A = 1/2.
c_V = 0.5
c_A = 0.5

# Step 1: Determine X2
# The problem gives the decay rate as:
# Gamma = (G_F * m_Z^3 / (12*sqrt(2)*pi)) * (c_V^2 + c_A^2)
# And defines X2 such that: Gamma = X2 * G_F * m_Z^3
# By comparing these, we find X2 = (c_V^2 + c_A^2) / (12*sqrt(2)*pi).

cv2_plus_ca2 = c_V**2 + c_A**2 # This is 0.5, or 1/2.

# Therefore, X2 = (1/2) / (12*sqrt(2)*pi) = 1 / (24*sqrt(2)*pi).
# We will use this symbolic form to maintain precision.
print(f"Step 1: Determine X2.")
print(f"Given c_V = {c_V} and c_A = {c_A}, the term (c_V^2 + c_A^2) is {cv2_plus_ca2}.")
print("From the provided formulas, we find X2 = (c_V^2 + c_A^2) / (12*sqrt(2)*pi).")
print("This simplifies to X2 = 1 / (24 * sqrt(2) * pi).\n")

# Step 2: Relate X1 and X2
# The general formula relating the decay rate (Gamma) and the spin-averaged
# squared amplitude (|M|^2) for a decay to two massless particles is:
# Gamma = |M|^2 / (16*pi*m_Z)
# Substituting the problem's definitions, |M|^2 = X1*G_F*m_Z^4 and Gamma = X2*G_F*m_Z^3:
# X2*G_F*m_Z^3 = (X1*G_F*m_Z^4) / (16*pi*m_Z)
# Simplifying this equation yields: X1 = 16*pi*X2.
print("Step 2: Relate X1 and X2.")
print("Using the kinematic relation Gamma = |M|^2 / (16*pi*m_Z), we find a simple relation:")
print("X1 = 16 * pi * X2\n")

# Step 3: Calculate (X1 * X2)^-1
# We substitute the expression for X1 into the quantity we want to calculate.
# (X1 * X2)^-1 = ((16*pi * X2) * X2)^-1 = (16*pi * X2^2)^-1
# Now we use the value of X2 we found:
# X2^2 = (1 / (24*sqrt(2)*pi))^2 = 1 / (24^2 * 2 * pi^2) = 1 / (1152 * pi^2)
# So, (16*pi * X2^2) = 16*pi / (1152*pi^2) = 16 / (1152*pi) = 1 / (72*pi)
# The inverse is 72*pi.
final_coefficient = 72
final_result_symbolic = f"{final_coefficient} * pi"
final_result_numeric = final_coefficient * math.pi

print("Step 3: Calculate the final result.")
print("(X1 * X2)^-1 = (16 * pi * X2^2)^-1")
print(f"= (16 * pi * [1 / (24*sqrt(2)*pi)]^2)^-1")
print(f"= (16 * pi / (1152 * pi^2))^-1")
print(f"= (1 / ({final_coefficient} * pi))^-1\n")

# The final equation as requested.
print(f"The final equation is (X1*X2)^-1 = {final_coefficient} * pi")
print(f"The numerical value is approximately: {final_result_numeric}")

print(f"<<<{final_result_symbolic}>>>")