from sympy import bernoulli, Rational

# Step 1: Define the formula for the Euler characteristic of the moduli stack M_3.
# chi(M_g) = -B_{2g} / (2g), where B_{2g} is the 2g-th Bernoulli number.
# For genus g = 3, we need B_6.
g = 3
B6 = bernoulli(2 * g)

# Calculate chi(M_3)
chi_M3 = -Rational(B6, 2 * g)

# Step 2: Use the known value for the Euler characteristic of the hyperelliptic locus H_3.
# From the literature (e.g., Bergström, Tommasi), chi(H_3) = -1/120.
chi_H3 = Rational(-1, 120)

# Step 3: The stack of smooth plane quartics [U/G] corresponds to the stack of
# non-hyperelliptic genus 3 curves, M_3^nh.
# The Euler characteristic is chi(M_3^nh) = chi(M_3) - chi(H_3).
result = chi_M3 - chi_H3

# Step 4: Print the final calculation step by step.
print("The problem is to compute the orbifold Euler characteristic of the quotient stack [U/G] of smooth plane quartics.")
print("This stack is the moduli stack of non-hyperelliptic genus 3 curves, M_3^nh.")
print("The Euler characteristic is computed using the decomposition of the moduli stack of genus 3 curves, M_3:")
print("chi([U/G]) = chi(M_3^nh) = chi(M_3) - chi(H_3)")
print("\nFirst, we calculate chi(M_3) using the Harer-Zagier formula chi(M_g) = -B_{2g}/(2g).")
print(f"For g=3, the 6th Bernoulli number B_6 is {B6}.")
print(f"So, chi(M_3) = -({B6}) / (2*3) = {chi_M3}.")
print("\nNext, we use the known value for the hyperelliptic locus, chi(H_3).")
print(f"chi(H_3) = {chi_H3}.")
print("\nFinally, we compute the difference:")
print(f"chi([U/G]) = chi(M_3) - chi(H_3) = ({chi_M3}) - ({chi_H3}) = {result}")
<<<11/2520>>>