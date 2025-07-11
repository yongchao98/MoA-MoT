import collections

def format_jones_polynomial(poly_dict):
    """
    Formats a dictionary of polynomial terms into a human-readable string.
    The dictionary keys are the powers and values are the coefficients.
    The output is ordered by decreasing power.
    """
    
    # Sort the polynomial terms by power in descending order
    sorted_terms = sorted(poly_dict.items(), key=lambda item: item[0], reverse=True)
    
    # Filter out terms with zero coefficients
    sorted_terms = [term for term in sorted_terms if term[1] != 0]

    if not sorted_terms:
        print("0")
        return

    result_parts = []
    
    # Process the first term separately to handle the leading sign
    first_power, first_coeff = sorted_terms[0]
    
    # Coefficient part
    coeff_str = ""
    if abs(first_coeff) != 1 or first_power == 0:
        coeff_str = str(abs(first_coeff))
    
    if first_coeff < 0:
        coeff_str = "-" + coeff_str

    # Variable part
    var_str = ""
    if first_power != 0:
        if first_power == 1:
            var_str = "t"
        else:
            var_str = f"t^{{{first_power}}}"

    result_parts.append(coeff_str + var_str)
    
    # Process the remaining terms
    for power, coeff in sorted_terms[1:]:
        # Sign part
        sign_str = " + " if coeff > 0 else " - "
        
        # Coefficient part
        coeff_str = ""
        if abs(coeff) != 1 or power == 0:
            coeff_str = str(abs(coeff))
        
        # Variable part
        var_str = ""
        if power != 0:
            if power == 1:
                var_str = "t"
            else:
                var_str = f"t^{{{power}}}"

        result_parts.append(sign_str + coeff_str + var_str)

    # Join all parts and format curly braces for exponents
    final_string = "".join(result_parts).replace("^{", "^").replace("}", "")
    print(final_string)

# The Jones polynomial for the 9_42 knot
# V(t) = t - 1 + t^{-1} - 2t^{-2} + 2t^{-3} - 2t^{-4} + 2t^{-5} - t^{-6} + t^{-7}
jones_poly = {
    1: 1,
    0: -1,
    -1: 1,
    -2: -2,
    -3: 2,
    -4: -2,
    -5: 2,
    -6: -1,
    -7: 1
}

format_jones_polynomial(jones_poly)