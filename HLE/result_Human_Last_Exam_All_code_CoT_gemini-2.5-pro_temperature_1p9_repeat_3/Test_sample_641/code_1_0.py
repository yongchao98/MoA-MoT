# A script to calculate the number of involutions in PSU(4, 997)

# Set the value of q
q = 997

# Calculate the size of the first class of involutions (N1)
# Formula: N1 = q^4 * (q^2 + 1) * (q^2 - q + 1)
n1_val = q**4 * (q**2 + 1) * (q**2 - q + 1)

# Calculate the size of the second class of involutions (N2)
# Formula: N2 = q^2 * (q^3 + 1)
n2_val = q**2 * (q**3 + 1)

# The total number of involutions is the sum of the sizes of the two classes
total_involutions = n1_val + n2_val

# As requested, output the numbers in the final equation.
# The final equation is Total = N1 + N2
print("The number of involutions is calculated by summing the sizes of the two conjugacy classes of involutions.")
print(f"Size of the first class (N1): {n1_val}")
print(f"Size of the second class (N2): {n2_val}")
print("\nFinal Equation:")
print(f"{n1_val} + {n2_val} = {total_involutions}")

# For clarity, let's also print the compact polynomial result.
# Total = q^8 - q^7 + 2*q^6 + q^4 + q^2
compact_total = q**8 - q**7 + 2*q**6 + q**4 + q**2
print(f"\nVerification using simplified polynomial: {compact_total}")