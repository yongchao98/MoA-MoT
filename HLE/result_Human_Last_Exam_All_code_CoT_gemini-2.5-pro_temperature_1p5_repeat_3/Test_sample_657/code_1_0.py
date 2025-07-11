print("The set of all proper stabilizing controllers H_2(s) for the plant H_1(s) = s/(s^2-1) is given by the Youla-Kucera parametrization:")
print("")
print("H_2(s) = N_c(s) / D_c(s)")
print("")
print("The numerator N_c(s) and denominator D_c(s) are parametrized by a function K(s), which must be any stable and proper rational function.")
print("Note: A constant value is a valid choice for K(s).")
print("")

# Define coefficients for the numerator expression
num_k_poly = "(1*s^2 + 0*s - 1)"
num_const_poly = "(4*s^2 + 8*s + 4)"
print(f"Numerator: N_c(s) = {num_k_poly} * K(s) + {num_const_poly}")

# Define coefficients for the denominator expression
den_k_poly = "(-1*s + 0)"
den_const_poly = "(1*s^2 + 0*s - 1)"
print(f"Denominator: D_c(s) = {den_k_poly} * K(s) + {den_const_poly}")
