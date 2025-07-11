def solve_scattering_statements():
    """
    Analyzes statements about multichannel quantum scattering.

    The core principle is that the nature of the coupling (diagonal, trivially
    coupled, or nontrivially coupled) is consistent across the potential V(r),
    the Jost matrix F(E), and the scattering matrix S(E).

    - If V(r) is trivially coupled, so are F(E) and S(E).
    - If S(E) is trivially coupled, inverse scattering uniqueness implies V(r) is
      trivially coupled, and thus F(E) is as well.

    This establishes an equivalence: V is nontrivial <=> F is nontrivial <=> S is nontrivial.

    Based on this, we evaluate the statements:
    1. Nontrivial S <=> Nontrivial V. (Correct)
    2. Diagonal S => Diagonal V. (Correct)
    3. Nontrivial V <=> Nontrivial F. (Correct)
    4. Nontrivial F <=> Nontrivial S. (Correct)
    5. Exists Nontrivial V with Diagonal F. (Incorrect, contradicts the equivalence).
    """
    correct_statements = [1, 2, 3, 4]
    
    # Create the equation string showing each correct statement number
    equation_str = " + ".join(map(str, correct_statements))
    
    # Calculate the sum
    total = sum(correct_statements)
    
    # Print the final equation as requested
    print(f"{equation_str} = {total}")

solve_scattering_statements()