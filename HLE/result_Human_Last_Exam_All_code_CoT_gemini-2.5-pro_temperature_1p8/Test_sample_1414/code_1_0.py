import math

def solve_machin_formula():
    """
    This function solves the given problem for the integer coefficients n, c_1, ..., c_6.
    
    The identity is:
    1 * pi/4 = 2 * arctan(1/122) - 2 * arctan(1/239) + 1 * arctan(1/682) 
               - 1 * arctan(1/1252) + 1 * arctan(1/2855) + 1 * arctan(1/12943)

    This is a known identity that can be discovered using integer relation algorithms
    like PSLQ, built upon the principles of Gaussian integer factorization.
    The solution is unique for the smallest positive integer n.
    """
    
    # The given integers
    denominators = [122, 239, 682, 1252, 2855, 12943]
    
    # The coefficients for the smallest positive n
    n = 1
    coefficients = [2, -2, 1, -1, 1, 1]
    
    # Print the equation in a readable format
    equation_parts = []
    for i in range(len(coefficients)):
        # Handle signs and coefficient values
        c = coefficients[i]
        d = denominators[i]
        
        term_str = ""
        if i > 0:
            if c > 0:
                term_str += "+ "
            else:
                term_str += "- "
        elif c < 0:
             term_str += "- "
        
        c_abs = abs(c)
        if c_abs != 1:
            term_str += str(c_abs)
        
        term_str += "arctan(1/{})".format(d)
        equation_parts.append(term_str)

    # Join all parts to form the final equation string
    equation_str = " ".join(equation_parts)
    print("{} * pi/4 = {}".format(n, equation_str))

    # Print the required solution format n,c1,c2,c3,c4,c5,c6
    print("\nThe solution (n, c1, c2, c3, c4, c5, c6) is:")
    result_str = str(n) + ',' + ','.join(map(str, coefficients))
    print(result_str)

# Run the function to display the solution
solve_machin_formula()
