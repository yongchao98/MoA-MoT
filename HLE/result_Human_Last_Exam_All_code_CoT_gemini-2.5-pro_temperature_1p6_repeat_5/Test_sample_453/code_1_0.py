from fractions import Fraction

def harmonic(k):
    """Calculates the k-th harmonic number as a Fraction for exact arithmetic."""
    if k <= 0:
        return Fraction(0)
    return sum(Fraction(1, i) for i in range(1, k + 1))

def calculate_energy_shift_coefficient(n, l):
    """
    Calculates the numerical coefficient for the second-order relativistic energy shift
    using the simplified formula for l = n - 1 states.
    The energy shift is of the form: ΔE = Coeff * m * c^2 * α^6.
    This function computes 'Coeff'.
    
    The formula is derived from results in advanced quantum mechanics literature 
    (e.g., Zapryagaev & Manakov, 1976), which simplifies for l = n-1 to:
    ΔE = (m*c^2*α^6 / n^3) * [ -1/(2*(l+1/2)^2) - S1 ]
    where S2 term is zero, and S1 is:
    S1 = 1/(n*(l+1/2)) - 2/n^2 + (2/(l+1/2) - 4/n) * (H_{n+l} - H_{l})
    (Using ψ(x+1) = ψ(x) + 1/x relation to H)
    A more careful conversion from digamma functions to harmonic numbers gives:
    ψ(2n) - ψ(n) = H_{2n-1} - H_{n-1}
    """
    if l != n - 1:
        raise ValueError("This simplified formula is only valid for l = n - 1.")

    # Convert to fractions for precision
    n_f = Fraction(n)
    l_f = Fraction(l)

    # Term S1 from the formula
    l_plus_half = l_f + Fraction(1, 2)
    s1_term1 = Fraction(1, n_f * l_plus_half)
    s1_term2 = -Fraction(2, n_f**2)
    
    # Harmonic number difference from digamma functions ψ(2n) - ψ(n) = H_{2n-1} - H_{n-1}
    h_diff = harmonic(2 * n - 1) - harmonic(n - 1)
    
    s1_h_coeff = Fraction(2, l_plus_half) - Fraction(4, n_f)
    s1_term3 = s1_h_coeff * h_diff

    s1 = s1_term1 + s1_term2 + s1_term3

    # Main term in the brackets
    main_term = -Fraction(1, 2 * l_plus_half**2)
    
    bracket_value = main_term - s1
    
    # Final coefficient
    coefficient = Fraction(1, n_f**3) * bracket_value
    
    return coefficient

# Given state parameters
n_val = 3
l_val = 2

# Calculate the coefficient
coeff = calculate_energy_shift_coefficient(n_val, l_val)
numerator = coeff.numerator
denominator = coeff.denominator

# Print the final result in the required format
print("The second-order energy shift ΔE for the n=3, l=2 state is given by the following equation:")
print(f"ΔE = ({numerator} / {denominator}) * m * c^2 * α^6")
print("\nWhere:")
print("m = mass of the electron")
print("c = speed of light")
print("α = fine-structure constant")
print(f"The final equation with the calculated numbers is:")
print(f"ΔE = ({numerator} / {denominator}) * m * c^2 * α^6")
