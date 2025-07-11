from fractions import Fraction

def solve_mass_problem(n, q):
    """
    Calculates the total mass based on the derived formula.

    The problem asks for the total mass of (q_v * (q-1) / (q_v-1)) * mu, where mu is the measure
    on the space GL_n^1(K) / GL_n(R).
    
    Based on the theory of volumes of S-arithmetic quotients, we can derive the formula for the mass.
    The setup strongly suggests K is the completion of a function field F=F_q(T) at a place,
    so the residue field size q_v is equal to q. The affine ring R is then F_q[T].

    The total mass of mu is given by:
    mass(mu) = (1 / (q-1)) * Product_{k=2 to n} Z_F(k)
    where Z_F(k) is the Dedekind zeta function of F = F_q(T) evaluated at k.
    Z_F(s) = 1 / ((1 - q^-s) * (1 - q^(1-s))).

    The factor in the question becomes q * (q-1) / (q-1) = q.
    So, the final quantity is: q * mass(mu) = (q / (q-1)) * Product_{k=2 to n} Z_F(k).
    
    Args:
        n (int): The dimension for GL_n.
        q (int): A prime power, the order of the residue field.
    
    Returns:
        The total mass as a Fraction object.
    """
    if not isinstance(n, int) or n < 1:
        raise ValueError("n must be a positive integer.")
    if not isinstance(q, int) or q < 2:
        raise ValueError("q must be a prime power >= 2.")

    # Prefactor for the final mass calculation
    prefactor = Fraction(q, q - 1)

    if n == 1:
        # The product is empty (equals 1)
        total_mass = prefactor
        print(f"For n=1, q={q}:")
        print(f"The space is O*/R* = O*/F_q^*.")
        print(f"The measure mu has mass 1/(q-1) = 1/{q-1}.")
        print(f"The requested quantity is q * mu = {q} * {Fraction(1, q-1)} = {total_mass}.")
        print(f"Final equation: {q} / ({q} - 1) = {total_mass}")
        return total_mass

    # Calculate zeta values and the total product
    zeta_values = []
    total_product = Fraction(1)
    
    for k in range(2, n + 1):
        # Z_F(k) = 1 / ((1 - q^-k) * (1 - q^-(k-1)))
        term1 = 1 - Fraction(1, q**k)
        term2 = 1 - Fraction(1, q**(k-1))
        zeta_k = 1 / (term1 * term2)
        zeta_values.append(zeta_k)
        total_product *= zeta_k

    total_mass = prefactor * total_product

    # Output the equation details
    print(f"Calculation for n={n}, q={q}:")
    equation_str = f"{q} / ({q} - 1)"
    for i, z in enumerate(zeta_values):
        k = i + 2
        print(f"Z_F({k}) = 1 / ((1 - {q}^-{k}) * (1 - {q}^-({k-1}))) = {z}")
        equation_str += f" * Z_F({k})"

    final_eq_str = f"{prefactor}"
    for z in zeta_values:
        final_eq_str += f" * {z}"
    
    print(f"\nFinal equation: {equation_str} = {total_mass}")
    print(f"Substituting values: {final_eq_str} = {total_mass}")
    
    return total_mass

if __name__ == '__main__':
    # Example values for n and q, as they are not specified in the problem.
    # The user can change these values.
    n_val = 4
    q_val = 3
    
    final_answer = solve_mass_problem(n_val, q_val)
    print(f"\nThe numerical result for n={n_val}, q={q_val} is {final_answer.numerator}/{final_answer.denominator}")
