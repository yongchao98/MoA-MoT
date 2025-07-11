import math

def calculate_sequences():
    """
    This function calculates the number of different score sequences possible.
    The problem is equivalent to finding the coefficient of x^48 in the expansion of:
    P(x) = (1+x+...+x^6)^8 * (1+x+...+x^8)^6
         = (1-x^7)^8 * (1-x^9)^6 * (1-x)^(-14)
    We use the binomial expansions of these terms to find the coefficient.
    """
    total_count = 0
    
    print("The final number is a sum of the following terms:")
    print("-" * 60)
    # Loop for the (1-x^7)^8 term, corresponding to the 8 players of team A
    for i in range(9):  # i from 0 to 8
        # Loop for the (1-x^9)^6 term, corresponding to the 6 players of team B
        for j in range(7):  # j from 0 to 6
            
            # We are looking for the coefficient of x^k in (1-x)^(-14)
            # such that 7*i + 9*j + k = 48
            k = 48 - 7 * i - 9 * j
            
            if k < 0:
                continue

            # Binomial coefficient for (1-x^7)^8 is C(8, i)
            # Binomial coefficient for (1-x^9)^6 is C(6, j)
            # Binomial coefficient for (1-x)^(-14) is C(k + 14 - 1, 14 - 1) = C(k + 13, 13)
            
            # The sign of the term is (-1)^i * (-1)^j = (-1)^(i+j)
            sign = (-1)**(i + j)
            
            c1 = math.comb(8, i)
            c2 = math.comb(6, j)
            c3 = math.comb(k + 13, 13)
            
            term = sign * c1 * c2 * c3
            
            if term == 0:
                continue

            # Output each number in the final equation
            sign_char = '+' if sign > 0 else '-'
            print(f"Term (i={i}, j={j}): {sign_char} C(8,{i})*C(6,{j})*C({k+13},13) = {sign_char} {c1}*{c2}*{c3} = {term}")
            
            total_count += term
            
    print("-" * 60)
    print(f"Total number of different sequences: {total_count}")

calculate_sequences()