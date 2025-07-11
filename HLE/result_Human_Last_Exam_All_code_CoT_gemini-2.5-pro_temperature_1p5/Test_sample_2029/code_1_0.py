import math

def solve():
    """
    This script calculates the number of different score sequences (a_1,...,a_8, b_1,...,b_6)
    in a tournament between two teams.

    The problem is equivalent to finding the number of integer solutions to the equation:
    a_1 + ... + a_8 + b_1 + ... + b_6 = 48
    with the constraints: 0 <= a_i <= 6 and 0 <= b_j <= 8.

    This is solved using generating functions. The answer is the coefficient of x^48 in the
    expansion of P(x) = (1+...+x^6)^8 * (1+...+x^8)^6.

    This coefficient can be calculated with the formula:
    Sum_{i=0 to 8, j=0 to 6} [ (-1)^(i+j) * C(8, i) * C(6, j) * C(61 - 7*i - 9*j, 13) ]
    where C(n, k) is the binomial coefficient "n choose k".
    """

    # A helper function for combinations that handles n < k cases gracefully.
    def combinations(n, k):
        if k < 0 or k > n:
            return 0
        return math.comb(n, k)

    total_sequences = 0
    equation_parts = []

    print("The total number of sequences is given by the sum of the following terms:")

    # Iterate through the terms of the expansion of (1-x^7)^8, indexed by i
    # and (1-x^9)^6, indexed by j.
    for i in range(9):
        for j in range(7):
            # Contribution to the power of x from the numerators (1-x^7)^8 and (1-x^9)^6
            power_from_numerator = 7 * i + 9 * j

            # We need the coefficient of x^48 in total. If this part already exceeds 48,
            # we don't need to consider it.
            if power_from_numerator > 48:
                continue
            
            # The coefficient from the (1-x)^-14 part needs to be for x^k where k = 48 - power_from_numerator.
            k = 48 - power_from_numerator
            
            # The full term in the sum is product of coefficients from the three parts:
            # (1-x^7)^8, (1-x^9)^6, and (1-x)^-14
            term_value = ((-1)**(i + j) * 
                          combinations(8, i) * 
                          combinations(6, j) * 
                          combinations(k + 13, 13))
            
            # Terms with a value of 0 don't need to be printed or added.
            if term_value == 0:
                continue

            # Add to the total count.
            total_sequences += term_value
            
            # Format the string for this part of the equation.
            sign = "+" if term_value > 0 else "-"
            # We show each number in the equation, as requested.
            # For C(n,k): n comes from the coefficient part and k from the summation index.
            # For the last C(n,k): n is k+13 and k is 13.
            equation_parts.append(f"{sign} C(8, {i}) * C(6, {j}) * C({k+13}, 13)")

    # Adjust the sign for the first term if it's positive.
    if equation_parts[0].startswith("+"):
        equation_parts[0] = equation_parts[0][2:]
    
    # Print the full equation.
    print("\n".join(equation_parts))
    
    # Print the final result.
    print(f"\nThe total number of different sequences is:\n{total_sequences}")

solve()