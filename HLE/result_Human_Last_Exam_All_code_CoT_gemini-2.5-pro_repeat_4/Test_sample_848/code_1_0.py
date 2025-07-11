import math

def solve():
    """
    Calculates the integer part of 10^4 * lim(F(N)/ln(N)) as N -> infinity.
    """
    # The two families of solutions are characterized by k=13 and k=5.
    k1 = 13
    k2 = 5

    # The growth rate of each family is determined by the larger root of
    # the characteristic equation lambda^2 - k*lambda + 1 = 0.
    
    # For k=13:
    d1 = k1**2 - 4
    sqrt_d1 = math.sqrt(d1)
    phi1 = (k1 + sqrt_d1) / 2
    log_phi1 = math.log(phi1)
    limit_term1 = 2 / log_phi1

    # For k=5:
    d2 = k2**2 - 4
    sqrt_d2 = math.sqrt(d2)
    phi2 = (k2 + sqrt_d2) / 2
    log_phi2 = math.log(phi2)
    limit_term2 = 2 / log_phi2
    
    # The final limit is the sum of the contributions from both families.
    L = limit_term1 + limit_term2
    
    # The problem asks for the integer part of 10^4 * L.
    final_value = 10000 * L
    
    # Outputting the numbers in the final equation as requested.
    print(f"The calculation is for the expression: 10^4 * (2 / ln((k1 + sqrt(d1))/2) + 2 / ln((k2 + sqrt(d2))/2))")
    print(f"For the first family (k1 = {k1}):")
    print(f"  d1 = {k1}^2 - 4 = {d1}")
    print(f"  sqrt(d1) = {sqrt_d1}")
    print(f"  phi1 = ({k1} + {sqrt_d1}) / 2 = {phi1}")
    print(f"  ln(phi1) = {log_phi1}")
    print(f"  Term 1 = 2 / {log_phi1} = {limit_term1}")
    print("")
    print(f"For the second family (k2 = {k2}):")
    print(f"  d2 = {k2}^2 - 4 = {d2}")
    print(f"  sqrt(d2) = {sqrt_d2}")
    print(f"  phi2 = ({k2} + {sqrt_d2}) / 2 = {phi2}")
    print(f"  ln(phi2) = {log_phi2}")
    print(f"  Term 2 = 2 / {log_phi2} = {limit_term2}")
    print("")
    print(f"The limit L = Term 1 + Term 2 = {L}")
    print(f"The final value is 10^4 * L = {final_value}")
    print(f"The integer part of the final value is: {math.floor(final_value)}")

solve()