import sympy

def calculate_proportionality_factor():
    """
    Calculates the proportionality factor C(d, k) for the gamma matrix identity.
    The user is prompted to enter the number of dimensions 'd' and the rank 'k'.
    """
    try:
        d_str = input("Enter the number of dimensions (d): ")
        k_str = input("Enter the rank of the gamma matrix product (k): ")
        
        d = sympy.sympify(d_str)
        k = sympy.sympify(k_str)

        if not (d.is_integer or d.is_symbol) or not (k.is_integer or k.is_symbol):
             raise ValueError("Inputs must be integers or symbolic variables.")
        if d.is_integer and k.is_integer:
            if k < 0 or d < 0:
                raise ValueError("Dimensions and rank must be non-negative.")
            if k > d:
                print(f"Warning: For k > d ({k} > {d}), the antisymmetrized gamma matrix is zero. The factor is irrelevant, but we calculate it formally.")
        
        # The proportionality factor C(d,k) based on literature results.
        # C(d, k) = 2 * (-1)**(k+1) * (d - 2*k)
        factor = 2 * ((-1)**(k + 1)) * (d - 2*k)
        
        print("\nThe gamma matrix product transforms as:")
        print(f"γ_μν * γ_{'μ₁...μₖ'} * γ^μν = C(d, k) * γ_{'μ₁...μₖ'}")
        
        print("\nThe proportionality factor C(d, k) is:")
        
        # To display the expression clearly
        term1 = 2
        term2 = (-1)**(k+1)
        term3 = (d-2*k)
        
        if term2.is_constant() and term2 == 1:
            if term1 == 1:
                 # Should not happen as term1 is 2
                 print(f"C({d}, {k}) = {term3}")
            else:
                 print(f"C({d}, {k}) = {term1} * ({term3}) = {factor}")
        elif term2.is_constant() and term2 == -1:
            print(f"C({d}, {k}) = -{term1} * ({term3}) = {factor}")
        else: # symbolic exponent
            print(f"C({d}, {k}) = {term1} * ({term2}) * ({term3}) = {factor}")

    except (ValueError, TypeError, sympy.SympifyError) as e:
        print(f"Invalid input: {e}")

if __name__ == "__main__":
    calculate_proportionality_factor()
