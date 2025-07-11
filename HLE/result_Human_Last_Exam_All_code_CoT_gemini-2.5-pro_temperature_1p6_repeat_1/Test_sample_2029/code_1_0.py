import math

def combinations(n, k):
    """
    Calculates the number of combinations C(n, k).
    Returns 0 if k < 0 or k > n.
    """
    if k < 0 or k > n:
        return 0
    return math.comb(n, k)

def solve_and_print():
    """
    Calculates the total number of different score sequences based on the
    derived combinatorial formula.
    """
    total_sequences = 0
    
    print("The total number of sequences is calculated by the sum S, where:")
    print("S = Sum_{i=0 to 8, j=0 to 6} (-1)^(i+j) * C(8,i) * C(6,j) * C(61 - 7i - 9j, 13)\n")
    print("Below are the non-zero terms in the sum, which constitute the final equation:")
    print("-" * 52)
    print("Term Expression                        Value")
    print("-" * 52)
    
    # Lists to store the parts of the final equation
    symbolic_terms = []
    numerical_terms = []

    # Iterate through all possible values for i and j
    for i in range(9):  # i from 0 to 8
        for j in range(7):  # j from 0 to 6
            k = 7 * i + 9 * j
            
            # A term is non-zero only if the combination is well-defined.
            if 61 - k >= 13:
                # Calculate components of the term
                c8i = combinations(8, i)
                c6j = combinations(6, j)
                c_rem = combinations(61 - k, 13)
                
                sign = (-1)**(i + j)
                term_value = sign * c8i * c6j * c_rem
                
                # A term might be zero if C(n,k) is zero, e.g., C(8,9)
                if term_value != 0:
                    total_sequences += term_value
                    
                    # Format the parts for printing the equation
                    sign_char = '+' if sign > 0 else '-'
                    symbolic_term = f"{sign_char} C(8,{i})*C(6,{j})*C({61-k},13)"
                    symbolic_terms.append(symbolic_term)
                    numerical_terms.append(term_value)
                    
                    # Print each number/term of the equation
                    print(f"{symbolic_term:<36} {term_value:+,}")

    print("-" * 52)
    print(f"The final result is the sum of these terms.")
    print(f"\nTotal number of sequences = {total_sequences:,}")

# Execute the calculation and print the results
solve_and_print()