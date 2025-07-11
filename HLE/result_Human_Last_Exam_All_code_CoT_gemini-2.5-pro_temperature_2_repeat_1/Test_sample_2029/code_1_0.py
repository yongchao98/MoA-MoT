import math

def nCr_safe(n, r):
    """
    A safe combination function C(n, r) that returns 0 if r < 0 or r > n.
    """
    if r < 0 or r > n:
        return 0
    return math.comb(n, r)

def solve_match_sequences():
    """
    Calculates the number of different score sequences possible in the tournament.
    """
    # Step 1: Calculate coefficients u_m for U(x) = (1 - x^7)^8 * (1 - x^9)^6
    # We only need coefficients up to m=48.
    # u[m] will store the coefficient of x^m.
    u = [0] * 49

    # Expand U(x) = (sum_i C(8,i)(-x^7)^i) * (sum_j C(6,j)(-x^9)^j)
    # The term for x^m is when 7i + 9j = m
    for i in range(9):  # i from 0 to 8
        for j in range(7):  # j from 0 to 6
            m = 7 * i + 9 * j
            if m <= 48:
                term = nCr_safe(8, i) * nCr_safe(6, j) * ((-1)**(i + j))
                u[m] += term

    # Step 2: Calculate the final sum: sum_{m=0 to 48} u_m * C(61 - m, 13)
    total_sequences = 0
    equation_parts = []
    numerical_terms = []

    for m in range(49):  # m from 0 to 48
        u_m = u[m]
        if u_m != 0:
            # Coefficient of x^(48-m) in (1-x)^(-14) is C(48-m+13, 13)
            v_n = nCr_safe(61 - m, 13)
            term = u_m * v_n
            total_sequences += term

            # Store parts for printing the full equation
            part_str = f"{u_m}*C({61 - m}, 13)"
            equation_parts.append(part_str)
            numerical_terms.append(term)
    
    print("The total number of sequences is the result of the following sum:")
    # Construct and print the equation in a readable format
    final_equation_str = equation_parts[0]
    for part in equation_parts[1:]:
        if part.startswith('-'):
            # For negative terms, use " - "
            final_equation_str += f" - {part[1:]}"
        else:
            # For positive terms, use " + "
            final_equation_str += f" + {part}"
    print(final_equation_str)
    
    print("\nWhich evaluates to:")
    final_numerical_str = str(numerical_terms[0])
    for num in numerical_terms[1:]:
         if num < 0:
            final_numerical_str += f" - {abs(num)}"
         else:
            final_numerical_str += f" + {num}"
    print(final_numerical_str)
    

    print(f"\nFinal Answer: {total_sequences}")

solve_match_sequences()
<<<35985858000629235946>>>