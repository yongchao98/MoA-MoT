import mpmath

# Set the precision for mpmath calculations. High precision is essential for
# the PSLQ algorithm to find the correct integer relation.
mpmath.mp.dps = 150

# The denominators of the arctan arguments from the problem statement
x_values = [122, 239, 682, 1252, 2855, 12943]

# We are trying to solve:
# c1*arctan(1/122) + ... + c6*arctan(1/12943) = n*pi/4
# This can be rewritten as an integer relation finding problem:
# c1*arctan(1/122) + ... + c6*arctan(1/12943) - n*(pi/4) = 0
# The vector of values for which we seek a relation is:
v = [mpmath.atan(mpmath.mpf(1) / x) for x in x_values]
v.append(mpmath.pi / 4)

# Use the PSLQ algorithm to find the integer coefficients [c1, ..., c6, -n].
# It returns a list of integers `coeffs` such that sum(coeffs[i]*v[i]) is nearly zero.
coeffs = mpmath.pslq(v)

# The last coefficient corresponds to pi/4. From our equation setup, this is -n.
# We need to find the smallest *positive* n. The PSLQ algorithm provides the smallest
# integer vector (in norm), so we just need to ensure the sign of n is positive.
# If the last coefficient is positive, our n would be negative. To fix this, we
# can multiply the entire relation by -1.
if coeffs[-1] > 0:
    coeffs = [-c for c in coeffs]

# Extract n and the c_i coefficients from the result
n = -coeffs[-1]
c_coeffs = coeffs[:-1]

# Print the final equation with the discovered coefficients.
print(f"The discovered relation is:")

# Print the left hand side: n * pi / 4
print(f"{n} * pi / 4 = ", end="")

# Print the right hand side, term by term, for nice formatting
is_first_term = True
for c, x in zip(c_coeffs, x_values):
    if c == 0:
        continue
    
    # Add a "+" or "-" sign for all but the first term
    if not is_first_term:
        if c > 0:
            print(" + ", end="")
        else:
            print(" - ", end="")
    # Handle a potential leading negative sign
    elif c < 0:
        print("-", end="")

    is_first_term = False
    
    c_abs = abs(c)

    # Print the coefficient if it's not 1, then the arctan term
    if c_abs != 1:
        print(f"{c_abs}*arctan(1/{x})", end="")
    else:
        print(f"arctan(1/{x})", end="")

print("\n") # Add a newline for separation

# The problem asks for the answer in a specific format.
answer_string = f"{n},{','.join(map(str, c_coeffs))}"
print(f"The constants are (n, c1, c2, c3, c4, c5, c6):")
print(answer_string)