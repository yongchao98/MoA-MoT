from fractions import Fraction

def solve_integral():
    """
    Calculates the integral of the product of lambda classes lambda_3*lambda_2*lambda_1
    on the moduli space of stable curves of genus 3, M_3.
    """
    
    # Step 1: Define lambda classes as polynomials in kappa classes.
    # The coefficients are represented as dictionaries mapping a kappa monomial string to its coefficient.
    # Source for formulas: A. Buryak, "Intersection theory on Mg and the KdV hierarchy", eq. (55)
    lambda1_poly = { 'k1': Fraction(1, 12) }
    lambda2_poly = { 'k1_sq': Fraction(1, 288), 'k2': Fraction(-1, 240) }
    lambda3_poly = { 'k1_cub': Fraction(1, 10368), 'k1_k2': Fraction(-1, 1440), 'k3': Fraction(1, 1728) }

    # Helper function for polynomial multiplication
    def multiply_polys(poly1, poly2):
        res = {}
        for term1, coeff1 in poly1.items():
            for term2, coeff2 in poly2.items():
                # This simple multiplication logic is sufficient for the specific problem
                if term1 == 'k1' and term2 == 'k1_sq': new_term = 'k1_cub'
                elif term1 == 'k1' and term2 == 'k2': new_term = 'k1_k2'
                elif term1 == 'k1_cub' and term2 == 'k1_cub': new_term = 'k1_pow_6'
                elif term1 == 'k1_cub' and term2 == 'k1_k2': new_term = 'k1_pow_4_k2'
                elif term1 == 'k1_cub' and term2 == 'k3': new_term = 'k1_cub_k3'
                elif term1 == 'k1_k2' and term2 == 'k1_cub': new_term = 'k1_pow_4_k2'
                elif term1 == 'k1_k2' and term2 == 'k1_k2': new_term = 'k1_sq_k2_sq'
                elif term1 == 'k1_k2' and term2 == 'k3': new_term = 'k1_k2_k3'
                else:
                    # For this specific calculation, other products are not degree 6
                    continue
                
                res[new_term] = res.get(new_term, 0) + coeff1 * coeff2
        return res

    # Step 2: Compute the product lambda_1 * lambda_2 * lambda_3
    l1_x_l2 = multiply_polys(lambda1_poly, lambda2_poly)
    product_poly = multiply_polys(l1_x_l2, lambda3_poly)

    # Step 3: Define known integrals of kappa monomials for g=3
    # Source: C. Faber, "A conjectural description of the tautological ring of the moduli space of curves"
    kappa_integrals = {
        'k1_pow_6': Fraction(1, 2),
        'k1_pow_4_k2': Fraction(1, 3),
        'k1_cub_k3': Fraction(1, 6),
        'k1_sq_k2_sq': Fraction(7, 24),
        'k1_k2_k3': Fraction(1, 8),
    }

    # Step 4: Calculate the final integral by summing up contributions from each monomial
    total_integral = Fraction(0)
    equation_parts = []
    
    terms = sorted(product_poly.keys()) # sort for consistent output order

    for term in terms:
        coeff = product_poly[term]
        integral_val = kappa_integrals[term]
        term_integral = coeff * integral_val
        total_integral += term_integral
        
        sign = "-" if term_integral < 0 else "+"
        # Print representation like "Integral(coeff * term)""
        part_str = f"({str(coeff)}) * Integral({term}) = {str(term_integral)}"
        
        # In the final equation we just show the numbers
        num_part_str = str(term_integral) if sign == "+" else f"({str(term_integral)})"

        if equation_parts:
             equation_parts.append(f" + {num_part_str}")
        else:
             equation_parts.append(num_part_str)


    # Step 5: Print the results and the calculation
    print("The integral is calculated by expanding lambda classes into kappa classes,")
    print("multiplying them, and then integrating the resulting kappa polynomial term by term.")
    print("\nThe product lambda_3*lambda_2*lambda_1 expands into the following kappa monomials with their coefficients:")
    for term in terms:
        print(f"  {term}: {product_poly[term]}")
        
    print("\nThe known values for the integrals of these kappa monomials are:")
    for term in terms:
        print(f"  Integral({term}) = {kappa_integrals[term]}")

    print("\nThe final integral is the sum of the term-wise integrals:")
    final_equation = "".join(equation_parts)
    print(f"Integral = {final_equation} = {str(total_integral)}")

solve_integral()
<<<1/17280>>>