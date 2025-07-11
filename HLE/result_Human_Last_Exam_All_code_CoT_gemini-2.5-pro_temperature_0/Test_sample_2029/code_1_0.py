import math

def combinations(n, k):
    """
    Calculates the binomial coefficient C(n, k).
    Returns 0 if k < 0 or k > n.
    """
    if k < 0 or k > n:
        return 0
    return math.comb(n, k)

def solve_match_sequences():
    """
    Calculates the number of different win sequences based on the
    generating function method.
    """
    total_count = 0
    equation_parts = []

    # The formula is: Sum_{i=0 to 8} Sum_{j=0 to 6} (-1)^(i+j) * C(8,i) * C(6,j) * C(61 - 7i - 9j, 13)
    for i in range(9):  # i from 0 to 8
        for j in range(7):  # j from 0 to 6
            
            n_comb3 = 61 - 7 * i - 9 * j
            k_comb3 = 13

            # Calculate the value of the term
            term_value = combinations(8, i) * combinations(6, j) * combinations(n_comb3, k_comb3)

            if term_value == 0:
                continue

            # Determine the sign and add to the total
            if (i + j) % 2 == 1:
                total_count -= term_value
                sign = "-"
            else:
                total_count += term_value
                sign = "+"

            # Build the string for the equation part
            part_str = f"C(8,{i})*C(6,{j})*C({n_comb3},{k_comb3})"
            equation_parts.append(f"{sign} {part_str}")

    # Format the final equation string for printing
    full_equation = " ".join(equation_parts)
    if full_equation.startswith("+ "):
        full_equation = full_equation[2:]

    print("The number of sequences is given by the following calculation:")
    print(f"{full_equation} = {total_count}")

solve_match_sequences()
<<<3358314105011264>>>