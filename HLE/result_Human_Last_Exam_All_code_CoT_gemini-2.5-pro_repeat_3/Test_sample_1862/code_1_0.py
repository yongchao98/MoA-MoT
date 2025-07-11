# Known values from the problem description for the first particle (full sphere)
d1_nm = 50
t1_s = 20

# Known values for the second particle (quarter sphere)
d2_nm = 250

# Step 1: Calculate the radii
r1_nm = d1_nm / 2
r2_nm = d2_nm / 2

# Step 2: The relationship for lithiation time is t_new = t_old * (r_new / (2 * r_old))
# This is derived from t ∝ V/A, where for a full sphere V/A = r/3
# and for a quarter sphere V/A = r/6.
# t_new / t_old = (V_new/A_new) / (V_old/A_old) = (r_new/6) / (r_old/3) = r_new / (2 * r_old)
t2_s = t1_s * (r2_nm / (2 * r1_nm))

# Step 3: Print the calculation and the final result.
print("Calculating the lithiation time for the quarter-sphere particle.")
print(f"The relationship is: t_new = t_old * (r_new / (2 * r_old))")
print("Plugging in the values:")
print(f"t_new = {t1_s} s * ({r2_nm} nm / (2 * {r1_nm} nm))")
print(f"t_new = {t1_s} * ({r2_nm / (2 * r1_nm)})")
print(f"t_new = {t2_s} s")
print("\nFinal Equation:")
print(f"{t1_s} * ({r2_nm} / (2 * {r1_nm})) = {t2_s}")

# Final Answer in the specified format
final_answer = t2_s
print(f'<<<{final_answer}>>>')