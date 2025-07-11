import cypari2
import sys

# Set a higher precision for the calculations
# Allocate a larger PARI stack to handle complex computations.
pari = cypari2.pari
pari.pari_init(2**30, 2**28) # (parisize, parismax)
pari.default("realprecision", 50)

# Step 1: Define the elliptic curve E from its Weierstrass equation
# y^2 + a1*x*y + a3*y = x^3 + a2*x^2 + a4*x + a6
# For y^2 + y = x^3 - x^2 - 10x - 20, the coefficients are:
# a1=0, a2=-1, a3=1, a4=-10, a6=-20
try:
    E = pari.ellinit([0, -1, 1, -10, -20])
except cypari2.errors.PariError as e:
    print(f"An error occurred initializing the elliptic curve: {e}")
    sys.exit(1)

# Step 2: Compute the rank r of the Mordell-Weil group E(Q)
# ellrank returns a vector; the first component is the rank.
try:
    r = pari.ellrank(E)[0]
except cypari2.errors.PariError as e:
    print(f"An error occurred computing the rank: {e}")
    sys.exit(1)
    
# Step 3: Define the cubic primitive Dirichlet characters of conductor 7
# The group of characters mod 7 has a generator of order 6.
# The cubic characters are the ones of order 3.
all_chars_mod_7 = pari.chargroup(7)
cubic_chars = [chi for chi in all_chars_mod_7 if pari.charorder(chi) == 3]
chi1 = cubic_chars[0]
chi2 = cubic_chars[1]

# Step 4: Compute the leading coefficients a and b of the twisted L-functions at s=1
# The lfun function computes the leading term of the Taylor series L^{(k)}(1)/k!
try:
    a = pari.lfun(E, 1, chi1)
    b = pari.lfun(E, 1, chi2)
except cypari2.errors.PariError as e:
    print(f"An error occurred computing the L-function values: {e}")
    sys.exit(1)

# Step 5: Calculate the final sum r + a + b
# Since chi2 is the complex conjugate of chi1, b will be the conjugate of a.
# Their sum a+b is therefore a real number.
total_sum = r + a + b

# Output the components of the equation as requested
print("The problem is to compute r + a + b.")
print(f"The rank of the elliptic curve is r = {r}")
print(f"The leading coefficient of L(E, s, chi1) at s=1 is a = {a}")
print(f"The leading coefficient of L(E, s, chi2) at s=1 is b = {b}")
print(f"\nThe final equation is:")
print(f"{r} + ({a}) + ({b}) = {total_sum}")

# Round the final result to four decimal places
# The result should be real, so we can convert to float for rounding.
final_answer = round(float(total_sum), 4)

print(f"\nThe sum r + a + b rounded to four decimal places is: {final_answer:.4f}")