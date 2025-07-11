import math

def solve_polynomial_roots():
    """
    This function finds and prints the four roots of the specified polynomial
    in increasing order.
    """
    # Based on the analysis, the four roots of the polynomial are:
    # sqrt(14), 2*sqrt(6), sqrt(34), and 2*sqrt(11).
    # We represent each root as a tuple (coefficient, number_under_radical)
    # to fulfill the requirement of outputting each number in the final equation form.
    root_definitions = [
        (1, 14),  # Represents sqrt(14)
        (2, 6),   # Represents 2*sqrt(6)
        (1, 34),  # Represents sqrt(34)
        (2, 11)   # Represents 2*sqrt(11)
    ]

    # Calculate the numerical value for each root to sort them.
    # The value is calculated as: coefficient * sqrt(number).
    roots = []
    for coeff, num in root_definitions:
        value = coeff * math.sqrt(num)
        roots.append({
            "coeff": coeff,
            "num": num,
            "value": value
        })

    # Sort the list of roots based on their numerical value.
    sorted_roots = sorted(roots, key=lambda r: r['value'])

    print("The four roots of the polynomial in increasing order are:")
    
    final_answer = []
    for root in sorted_roots:
        c = root['coeff']
        n = root['num']
        v = root['value']
        final_answer.append(v)
        
        # The 'final equation' for each root is X = expression.
        # We print the numbers that form this expression.
        if c == 1:
            print(f"X = sqrt({n})  (approx. {v})")
        else:
            print(f"X = {c}*sqrt({n})  (approx. {v})")

    # The final answer in the specified format.
    print("\n<<<[3.7416573867739413, 4.898979485566356, 5.830951894845301, 6.6332495807108]>>>")

solve_polynomial_roots()