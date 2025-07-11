import cypari2
import math

# Initialize the cypari2 interface to the PARI/GP calculator
pari = cypari2.Pari()

# Step 1: Define the elliptic curve E.
# The user's curve is isogenous to '49a1' and they share the same L-function.
# We use the Cremona label '49a1' which is the standard minimal model.
E = pari.ellinit("49a1")

# Step 2: Determine the rank r of E(Q).
# The rank of the curve '49a1' is known to be 0.
r = 0

# Step 3: Define the cubic primitive Dirichlet characters of conductor 7.
# The character group of (Z/7Z)^* is cyclic of order 6. The cubic characters are the two elements of order 3.
# We find a generator of the character group and take its 2nd and 4th powers.
G = pari.chargroup(7)
gen = None
for char in G:
    if pari.zncharorder(char) == 6:
        gen = char
        break
chi1 = gen**2
chi2 = gen**4

# Step 4 & 5: Compute the leading coefficients a and b.
# The root number for both twists is -1, which implies the order of vanishing
# of the L-function at s=1 is 1. Thus, the leading coefficient is the
# first derivative L'(E, 1, chi).
# 'a' is the leading coefficient for L(E, s, chi1)
# 'b' is the leading coefficient for L(E, s, chi2)
a = pari.lfun(E, chi1, 1, 1)
b = pari.lfun(E, chi2, 1, 1)

# Step 6: Calculate the final sum r + a + b.
# Since b is the complex conjugate of a, their sum is 2 * Re(a).
final_value = r + a + b

# Print the final equation with all its components as requested.
# We display the complex numbers a and b and the final real result.
# The cypari complex number has real and imag parts which we can extract.
# We format the numbers for clarity in the equation.
a_str = f"({a.real()} + {a.imag()}*I)"
b_str = f"({b.real()} + {b.imag()}*I)"
final_value_str = f"{final_value.real()}"

print("The final equation is:")
print(f"{r} + {a_str} + {b_str} = {final_value_str}")

# Round the final result to four decimal places for the final answer.
rounded_result = round(final_value.real(), 4)

print(f"\nThe value of r + a + b rounded to four decimal places is: {rounded_result}")

# The final answer in the required format.
# print(f"<<<{rounded_result}>>>")