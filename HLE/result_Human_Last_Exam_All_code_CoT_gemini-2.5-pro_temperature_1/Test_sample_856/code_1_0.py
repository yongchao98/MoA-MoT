import math

def get_prime_factorization(n):
    factors = {}
    d = 2
    temp = n
    while d * d <= temp:
        while temp % d == 0:
            factors[d] = factors.get(d, 0) + 1
            temp //= d
        d += 1
    if temp > 1:
        factors[temp] = factors.get(temp, 0) + 1
    return factors

def phi(n, factors):
    result = n
    for p in factors:
        result = result * (p - 1) // p
    return result

def num_solutions_x2_1(factors):
    """Calculates the number of solutions to x^2 = 1 (mod n)"""
    num = 1
    for p, a in factors.items():
        if p == 2:
            if a == 1:
                num *= 1
            elif a == 2:
                num *= 2
            else: # a >= 3
                num *= 4
        else: # odd prime
            num *= 2
    return num

def num_solutions_x2_neg_1(factors):
    """Calculates the number of solutions to x^2 = -1 (mod n)"""
    num = 1
    for p, a in factors.items():
        if p == 2:
            if a > 1:
                return 0
        elif (p - 1) % 4 != 0:
            return 0
        else: # (p-1)%4 == 0
            num *= 2
    return num

def main():
    n = math.factorial(10)
    factors = get_prime_factorization(n)

    phi_n = phi(n, factors)
    num_x2_1 = num_solutions_x2_1(factors)
    num_x2_neg_1 = num_solutions_x2_neg_1(factors)
    
    # Number of solutions to 2x=0 in the group of units is 0 since n is even.
    num_2x_0 = 0
    
    # Using Burnside's Lemma for the group action q -> +/- q, q -> +/- q^{-1}
    # which is isomorphic to V_4 = C_2 x C_2
    # The number of orbits is 1/4 * (|U| + |{q:q=-q}| + |{q:q=q^-1}| + |{q:q=-q^-1}|)
    num_manifolds = (phi_n + num_2x_0 + num_x2_1 + num_x2_neg_1) // 4
    
    print(f"For the cyclic group C(10!), n = 10! = {n}")
    print(f"The number of units is phi(10!) = {phi_n}")
    print(f"Number of solutions to q^2 = 1 (mod 10!) is {num_x2_1}")
    print(f"Number of solutions to q^2 = -1 (mod 10!) is {num_x2_neg_1}")
    print(f"Number of solutions to 2q = 0 (mod 10!) in the group of units is {num_2x_0}")
    print("\nThe number of lens spaces is given by (phi(n) + num_sols(q^2=1) + num_sols(q^2=-1) + num_sols(2q=0)) / 4")
    print(f"Number of manifolds = ({phi_n} + {num_x2_1} + {num_x2_neg_1} + {num_2x_0}) / 4 = {num_manifolds}")

if __name__ == "__main__":
    main()
