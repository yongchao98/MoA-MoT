import math

def combinations(n, k):
    """
    Custom combinations function that returns 0 if k > n or k < 0.
    math.comb would raise an error.
    """
    if k < 0 or k > n:
        return 0
    return math.comb(n, k)

def solve_match_sequences():
    """
    Calculates the number of different score sequences (a_1,...,a_8,b_1,...,b_6)
    using the generating function method as outlined in the plan.
    """
    total_count = 0
    equation_terms = []

    # i is the power index from the (1-z^7)^8 expansion
    for i in range(9):  # C(8,i) is 0 for i > 8
        # j is the power index from the (1-z^9)^6 expansion
        for j in range(7):  # C(6,j) is 0 for j > 6
            
            # The power of z from the first two polynomial parts is 7i + 9j
            power_sum = 7 * i + 9 * j

            if power_sum <= 48:
                # We need the coefficient of z^(48 - power_sum) from (1-z)^(-14)
                n = 48 - power_sum
                
                # The coefficient is C(n + 13, 13)
                comb_c = combinations(n + 13, 13)
                
                # Coefficients from the first two polynomials
                comb_a = combinations(8, i)
                comb_b = combinations(6, j)
                
                term_value = comb_a * comb_b * comb_c
                
                # If term is zero, it contributes nothing and can be skipped
                if term_value == 0:
                    continue

                # The sign comes from (-1)^i and (-1)^j in the expansions
                sign = (-1)**(i + j)
                
                # Add the signed value to the total count
                total_count += sign * term_value
                
                # Store the signed value to print the full equation later
                equation_terms.append(sign * term_value)
                
    # As requested, print the equation showing each number being summed
    print("The total number of sequences is the result of the following sum:")
    
    # Format the equation string: term_1 +/- term_2 +/- ...
    equation_str = ""
    for k, term in enumerate(equation_terms):
        if k == 0:
            equation_str += str(term)
        else:
            if term >= 0:
                equation_str += f" + {term}"
            else:
                # For negative terms, print " - " followed by the positive value
                equation_str += f" - {-term}"
                
    print(f"{equation_str} = {total_count}")

# Execute the function
solve_match_sequences()