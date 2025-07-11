def solve_a4():
    """
    Calculates a(4), the maximal number of prime implicants of a Boolean function of 4 variables.

    This is achieved by the 4-variable parity function, whose prime implicants are its minterms.
    This script generates those minterms, forms the Boolean equation, and counts them.
    """
    n = 4
    variables = ['w', 'x', 'y', 'z']
    prime_implicants = []

    # Iterate through all 2^n possible minterms for n variables.
    # A minterm is represented by an integer from 0 to 2^n - 1.
    for i in range(2**n):
        # The parity function consists of minterms with an odd number of uncomplemented variables.
        # This corresponds to integers with an odd number of '1's in their binary representation.
        if bin(i).count('1') % 2 == 1:
            # This minterm is a prime implicant for the parity function.
            # We format it into a human-readable string.
            term_parts = []
            binary_representation = format(i, f'0{n}b')
            for j in range(n):
                if binary_representation[j] == '1':
                    term_parts.append(variables[j])  # Uncomplemented variable
                else:
                    term_parts.append(variables[j] + "'")  # Complemented variable
            prime_implicants.append("".join(term_parts))

    print("The Boolean function with the maximal number of prime implicants for n=4 is the parity function.")
    print("Its prime implicants are the minterms that make the function true.")
    
    # The final equation is the sum (OR) of these prime implicants.
    # We print each term of the equation as requested.
    equation = " + ".join(prime_implicants)
    print("\nThe final equation is: f(w,x,y,z) =", equation)

    # The value a(4) is the number of terms in this equation.
    result = len(prime_implicants)
    print(f"\nThe number of terms (prime implicants) in the equation is {result}.")
    print(f"Therefore, a(4) = {result}.")

solve_a4()