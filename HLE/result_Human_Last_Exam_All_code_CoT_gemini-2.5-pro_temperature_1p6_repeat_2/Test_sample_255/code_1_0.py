import math

def euler_totient(n):
    """
    Computes Euler's totient function phi(n).
    This counts the number of positive integers up to a given integer n
    that are relatively prime to n.
    """
    if n == 1:
        return 1
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
    Solves the problem by calculating the dimension of the cohomology group H^2(G,M).
    The dimension is given by the formula: phi(2) + phi(4) + phi(8).
    """
    print("The dimension of the cohomology group H^2(G,M) is given by the sum of Euler's totient function values.")
    print("Dimension = phi(2) + phi(4) + phi(8)")
    
    val2 = euler_totient(2)
    val4 = euler_totient(4)
    val8 = euler_totient(8)
    
    total_dimension = val2 + val4 + val8
    
    print(f"phi(2) = {val2}")
    print(f"phi(4) = {val4}")
    print(f"phi(8) = {val8}")
    print("So, the final dimension is:")
    print(f"{val2} + {val4} + {val8} = {total_dimension}")
    return total_dimension

if __name__ == "__main__":
    solve_cohomology_dimension()
