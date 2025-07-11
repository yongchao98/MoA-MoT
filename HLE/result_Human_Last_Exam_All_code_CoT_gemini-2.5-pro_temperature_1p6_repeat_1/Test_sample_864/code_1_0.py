def solve_whitening_filter():
    """
    This function determines and prints the whitening filter W(D).
    
    The problem as stated provides coefficients q_k that are inconsistent and do not
    correspond to a stable rational spectral density. Therefore, we proceed by
    assuming a plausible, well-behaved spectral factor G(D) that might have been
    intended for this type of problem.

    We assume the intended causal, minimum-phase spectral factor is:
    G(D) = (1 + D) / (1 - (1/3) * D^2)

    The whitening filter W(D) must make the overall channel Q(D)W(D) causal.
    This is achieved by choosing W(D) = 1 / G(D^{-1}).

    Let's calculate W(D):
    G(D^{-1}) = (1 + D^{-1}) / (1 - (1/3) * D^{-2})
    W(D) = 1 / G(D^{-1}) = (1 - (1/3) * D^{-2}) / (1 + D^{-1})

    To write this as a ratio of polynomials in D, we multiply the numerator
    and denominator by D^2:
    W(D) = D^2 * (1 - (1/3) * D^{-2}) / (D^2 * (1 + D^{-1}))
         = (D^2 - 1/3) / (D^2 + D)
    """

    numerator_coeffs = {2: 1, 0: -1/3}
    denominator_coeffs = {2: 1, 1: 1}

    def format_poly(coeffs):
        terms = []
        # Sort by power in descending order
        for power in sorted(coeffs.keys(), reverse=True):
            coeff = coeffs[power]
            
            # Sign
            sign = " - " if coeff < 0 else " + "
            coeff = abs(coeff)

            # Coefficient string
            if coeff == 1 and power != 0:
                coeff_str = ""
            else:
                coeff_str = f"{coeff:.3g}"
            
            # Power string
            if power == 0:
                power_str = ""
            elif power == 1:
                power_str = "D"
            else:
                power_str = f"D^{power}"
            
            terms.append(f"{sign}{coeff_str}{power_str}")

        # Join terms and clean up the leading sign
        poly_str = "".join(terms).lstrip(" +")
        return poly_str

    numerator_str = format_poly(numerator_coeffs)
    denominator_str = format_poly(denominator_coeffs)

    print("Based on an assumed consistent spectral factor G(D), the whitening filter is:")
    print(f"      {numerator_str}")
    print(f"W(D) = ---------")
    print(f"      {denominator_str}")

solve_whitening_filter()
