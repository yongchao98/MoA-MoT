import math

def solve_hamiltonian_period():
    """
    Solves the problem by finding n1, n2, and then calculating the period of the specified Hamiltonian.
    """
    # Step 1: Find n1 and n2.
    # The condition u_r(n) = n/2 - 1, combined with u_r(n +/- 1) = n +/- 1 - 1,
    # is satisfied by even integers n >= 4. This is because for even n, the potential
    # V(q) is a function of q^2, leading to an order reduction from n-1 to n/2 - 1.
    # For odd degrees (n+/-1), the standard order n-1 applies.
    # The 1st and 2nd smallest positive integers that fit this are 4 and 6.
    n1 = 4
    n2 = 6

    # Step 2: Calculate the argument of T, which is alpha.
    alpha_numerator = n1 - 1
    alpha_denominator = n2
    
    # Step 3: Calculate the period of the Hamiltonian H(p, q).
    # The Hamiltonian is given by H = 1/2 * (p^2 + q^2 - C * q^(n1/2)), where
    # C = (2/n1) * sqrt((n2 - n1) / (n1 / 2)).
    
    # Calculate the coefficient C
    C = (2 / n1) * math.sqrt((n2 - n1) / (n1 / 2))
    
    # The potential energy V(q) is the part of the Hamiltonian that depends on q.
    # From H = p^2/2 + V(q), we have:
    # V(q) = 1/2 * (q^2 - C * q^(n1/2))
    # Substitute n1 = 4 and the calculated value of C:
    # V(q) = 1/2 * q^2 - (C / 2) * q^2 = (1/2 - C / 2) * q^2
    
    # The potential is of the form V(q) = (1/2) * k * q^2 for a simple harmonic oscillator.
    # By comparing forms, k = 2 * (1/2 - C / 2) = 1 - C.
    spring_constant_k = 1 - C

    # The period T of a simple harmonic oscillator with Hamiltonian H = p^2/(2m) + (1/2)kq^2
    # is T = 2*pi*sqrt(m/k). In this problem, the mass m=1.
    # T = 2 * pi / sqrt(k)
    period = 2 * math.pi / math.sqrt(spring_constant_k)

    # Step 4: Output the final equation with all numbers.
    # The equation is T((n1-1)/n2) = period.
    print(f"T(({alpha_numerator}) / {alpha_denominator}) = {period}")

solve_hamiltonian_period()