import math

def solve():
    """
    Calculates the number of different score sequences (a_1,...,a_8, b_1,...,b_6).

    A score sequence is possible if and only if the sum of wins for Team A plus the sum of wins for Team B equals the total number of games, which is 48.
    Let S_A = a_1 + ... + a_8 and S_B = b_1 + ... + b_6. The constraint is S_A + S_B = 48.
    The number of sequences is the sum over all possible total scores k for Team A (from 0 to 48) of:
    (Number of ways Team A can score k) * (Number of ways Team B can score 48-k).

    This can be calculated using generating functions. The number is the coefficient of x^48 in the expansion of:
    F(x) = (1 + x + ... + x^6)^8 * (1 + x + ... + x^8)^6
    This is not quite right.
    
    The number of sequences is the constant term of P_A(x) * P_B(1/x) where
    P_A(x) = (1 + ... + x^6)^8 and P_B(x) = (1 + ... + x^8)^6.
    This leads to finding the coefficient of x^48 in the expansion of:
    (1-x^7)^8 * (1-x^9)^6 * (1-x)^(-14).

    The coefficient is given by the sum:
    Sum_{i,j} (-1)^(i+j) * C(8,i) * C(6,j) * C(48 - 7i - 9j + 13, 13)
    over valid integer ranges for i and j.
    """

    print("The total number of sequences is calculated by the following summation, where C(n,k) is the binomial coefficient 'n choose k':")

    total_sum = 0
    equation_parts = []

    # Iterate over j for the (1-x^9)^6 expansion. Max j is 5 since 9*6 > 48.
    for j in range(6):
        # Iterate over i for the (1-x^7)^8 expansion.
        i_limit = (48 - 9 * j) // 7
        for i in range(i_limit + 1):
            # Remaining power needed for the (1-x)^(-14) expansion.
            k = 48 - 7 * i - 9 * j
            
            # Calculate binomial coefficients for the term.
            comb_8_i = math.comb(8, i)
            comb_6_j = math.comb(6, j)
            # math.comb(n, k) is 0 if k > n.
            comb_k_13 = math.comb(k + 13, 13)
            
            # Calculate the value of the term.
            term_value = ((-1)**(i + j)) * comb_8_i * comb_6_j * comb_k_13
            total_sum += term_value
            
            # Format the string for the equation.
            if term_value >= 0:
                sign_str = "+ "
            else:
                sign_str = "- "

            if len(equation_parts) == 0 and sign_str == "+ ":
                sign_str = ""

            # Build the string showing each number in the term.
            part = f"{sign_str}{comb_8_i} * {comb_6_j} * {comb_k_13}"
            equation_parts.append(part)

    # Print the full equation.
    print(' '.join(equation_parts))
    
    # Print the final result.
    print(f"\n= {total_sum}")

solve()