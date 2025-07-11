import math

def main():
    """
    This program calculates the number of different sequences of win counts
    (a_1,...,a_8, b_1,...,b_6) that can occur in a match between
    Team A (8 members) and Team B (6 members).

    The problem is equivalent to finding the number of integer solutions to:
    a_1 + ... + a_8 + b_1 + ... + b_6 = 48
    subject to the constraints:
    0 <= a_i <= 6 (each member of Team A plays 6 games)
    0 <= b_j <= 8 (each member of Team B plays 8 games)

    This is solved using generating functions. The number is the coefficient
    of x^48 in the expansion of P(x) = (1+...+x^6)^8 * (1+...+x^8)^6.

    This simplifies to finding the coefficient of x^48 in:
    P(x) = (1-x^7)^8 * (1-x^9)^6 * (1-x)^(-14)

    The final formula for the coefficient is:
    Sum_{i=0 to 8, j=0 to 6} [ C(8,i)*(-1)^i * C(6,j)*(-1)^j * C(48-7i-9j+13, 13) ]
    """

    def combinations(n, k):
        # A robust combinations function that returns 0 if k > n or k < 0.
        if k < 0 or k > n:
            return 0
        return math.comb(n, k)

    total_sum = 0
    equation_parts = []

    # Iterate through the expansion of (1-x^7)^8 and (1-x^9)^6
    # i corresponds to the term in the expansion of (1-x^7)^8
    # j corresponds to the term in the expansion of (1-x^9)^6
    for i in range(8 + 1):
        for j in range(6 + 1):
            # We are looking for the coefficient of x^k in (1-x)^(-14)
            # where k = 48 - 7*i - 9*j
            k = 48 - (7 * i) - (9 * j)

            if k >= 0:
                # Coefficient from (1-x^7)^8 expansion term
                coeff_a = combinations(8, i)

                # Coefficient from (1-x^9)^6 expansion term
                coeff_b = combinations(6, j)
                
                # Sign from (-1)^i and (-1)^j
                sign = (-1)**(i + j)

                # Coefficient from (1-x)^(-14) expansion for x^k
                # This is C(k + 14 - 1, 14 - 1) = C(k + 13, 13)
                coeff_c = combinations(k + 13, 13)
                
                term = sign * coeff_a * coeff_b * coeff_c

                if term != 0:
                    total_sum += term
                    
                    # This part generates a string for the equation as requested
                    term_sign_str = " + " if term > 0 else " - "
                    
                    # Combine coefficients for a cleaner look in the equation
                    term_coeff_val = coeff_a * coeff_b
                    
                    term_str = f"{term_coeff_val}*C({k+13}, 13)"
                    
                    # Handle the first term of the equation
                    if not equation_parts:
                         # Remove leading " + " for the first term
                        if term_sign_str == " + ":
                            equation_parts.append(term_str)
                        else:
                            equation_parts.append(f"-{term_str}")
                    else:
                        equation_parts.append(term_sign_str + term_str)

    # To keep the output readable, we'll show a summary of the equation
    # and not the full, very long string.
    print("The calculation is based on the formula:")
    print("Sum_{i=0..8, j=0..6} [ C(8,i)*C(6,j)*(-1)^(i+j) * C(61-7i-9j, 13) ]")
    print("\nThe first few terms of the calculation are:")
    print(f"{' '.join(equation_parts[:4])} ...")
    
    print("\nThe final number of different sequences is:")
    print(total_sum)

if __name__ == "__main__":
    main()
