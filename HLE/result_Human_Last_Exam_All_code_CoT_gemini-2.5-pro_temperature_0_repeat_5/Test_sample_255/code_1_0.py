import math

def phi(n):
    """
    Computes Euler's totient function phi(n), which counts the number of positive integers up to a given integer n
    that are relatively prime to n.
    """
    if not isinstance(n, int) or n < 1:
        raise ValueError("Input must be a positive integer.")
    
    result = n
    p = 2
    temp_n = n
    while p * p <= temp_n:
        if temp_n % p == 0:
            while temp_n % p == 0:
                temp_n //= p
            result -= result // p
        p += 1
    if temp_n > 1:
        result -= result // temp_n
    return result

def solve_cohomology_dimension():
    """
    Calculates the dimension of the cohomology group H^2(G,M) based on the derived formula.
    The dimension is the sum of phi(d) for d=2, 4, 8.
    """
    divisors = [2, 4, 8]
    
    print("The dimension of the cohomology group H^2(G,M) is the sum of Euler's totient function phi(d) for d = 2, 4, 8.")
    print("-" * 30)

    dimensions = []
    for d in divisors:
        phi_d = phi(d)
        dimensions.append(phi_d)
        print(f"phi({d}) = {phi_d}")

    total_dimension = sum(dimensions)
    
    equation_parts = [str(d) for d in dimensions]
    equation = " + ".join(equation_parts)

    print("-" * 30)
    print(f"The dimension is the sum of these values:")
    print(f"Dimension = {equation} = {total_dimension}")
    
    return total_dimension

if __name__ == '__main__':
    solve_cohomology_dimension()