def print_chromatic_polynomial():
    """
    This function prints the final chromatic polynomial for the given graph.
    The polynomial P(k) = k^5 - 7k^4 + 19k^3 - 23k^2 + 10k was derived
    using the deletion-contraction method. This script formats the
    equation to explicitly show each coefficient.
    """
    # The coefficients for the polynomial, from the highest power (k^5) down to k^1.
    # The constant term (for k^0) is 0.
    coeffs = [1, -7, 19, -23, 10]
    
    terms = []
    # Iterate through the coefficients to build each term of the polynomial string.
    for i, c in enumerate(coeffs):
        power = 5 - i
        
        # Determine the sign and coefficient string for the term.
        if i == 0:
            # First term doesn't need a preceding sign unless it's negative.
            if c < 0:
                sign = "-"
                coeff_str = str(abs(c))
            else:
                sign = ""
                coeff_str = str(c)
        else:
            # Subsequent terms get a '+' or '-' sign.
            if c < 0:
                sign = " - "
                coeff_str = str(abs(c))
            else:
                sign = " + "
                coeff_str = str(c)

        # Determine the variable part of the term (e.g., k^5, k^4, k).
        if power > 1:
            var_str = f"k^{power}"
        else: # power is 1
            var_str = "k"
            
        # Assemble the full term string.
        terms.append(f"{sign}{coeff_str}*{var_str}")
        
    # Join all the terms to form the final equation string.
    equation = "".join(terms)
    
    print("The chromatic polynomial of the graph is P(k), where:")
    print(f"P(k) = {equation}")

# Execute the function to print the result.
print_chromatic_polynomial()