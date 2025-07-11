import mpmath

def solve_arctan_relation():
    """
    Finds the integer coefficients for the given Machin-like formula for pi.
    """
    # Set a high precision for the numerical calculations.
    # This is crucial for the PSLQ algorithm to find the correct integer relation.
    mpmath.mp.dps = 200

    # These are the denominators in the arctan terms from the problem statement.
    x_values = [122, 239, 682, 1252, 2855, 12943]

    # We construct a vector of the real numbers whose relation we want to find.
    # The vector is [pi/4, arctan(1/122), arctan(1/239), ...]
    reals_vector = [mpmath.pi / 4]
    for x in x_values:
        reals_vector.append(mpmath.atan(mpmath.mpf(1) / x))

    # The pslq function finds a non-zero integer vector 'coeffs' = [k0, k1, ..., k6]
    # such that k0*reals_vector[0] + k1*reals_vector[1] + ... = 0 (to a high precision).
    # Comparing with our equation n*(pi/4) - c1*arctan(...) - ... = 0, we have:
    # k0 = n
    # k1 = -c1, k2 = -c2, ..., k6 = -c6
    relation_coeffs = mpmath.pslq(reals_vector)

    # Extract n and the c_i coefficients from the result of PSLQ.
    n = relation_coeffs[0]
    c = [-k for k in relation_coeffs[1:]]

    # The problem requires the smallest *positive* n.
    # If PSLQ returns a relation where n is negative, we can multiply
    # the entire relation by -1 to get an equivalent valid relation.
    if n < 0:
        n = -n
        c = [-val for val in c]

    # Now we print the final result as requested.
    print("The solution is given by the following equation:")

    # Build and print the equation string, showing all the numbers involved.
    rhs_parts = []
    for i in range(len(c)):
        if c[i] == 0:
            continue

        # For the first term in the sum
        if not rhs_parts:
            # Handle coefficients of 1, -1, and others
            if c[i] == 1:
                rhs_parts.append(f"arctan(1/{x_values[i]})")
            elif c[i] == -1:
                rhs_parts.append(f"-arctan(1/{x_values[i]})")
            else:
                rhs_parts.append(f"{c[i]}*arctan(1/{x_values[i]})")
        # For subsequent terms
        else:
            sign = " + " if c[i] > 0 else " - "
            coeff_val = abs(c[i])
            # Handle coefficients of 1 and others
            if coeff_val == 1:
                rhs_parts.append(f"{sign}arctan(1/{x_values[i]})")
            else:
                rhs_parts.append(f"{sign}{coeff_val}*arctan(1/{x_values[i]})")

    rhs_string = "".join(rhs_parts)
    print(f"{n} * pi / 4 = {rhs_string}")

    # Print the final constants as a list
    print("\nThe values of the constants are:")
    result_string = f"{n},{c[0]},{c[1]},{c[2]},{c[3]},{c[4]},{c[5]}"
    print(result_string)

if __name__ == '__main__':
    solve_arctan_relation()