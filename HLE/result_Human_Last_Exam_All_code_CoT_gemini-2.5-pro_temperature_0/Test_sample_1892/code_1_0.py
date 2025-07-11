import math

def gcd(a, b):
    """Computes the greatest common divisor of a and b."""
    return math.gcd(a, b)

def phi(n):
    """Computes Euler's totient function."""
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

def solve_problem():
    """
    Calculates the parameters alpha and beta for the asymptotic formula.
    """
    k = 12
    
    # Step 1: Find phi(k)
    phi_k = phi(k)
    
    # Step 2: Find the residue classes a mod k where gcd(a, k) = 1
    residues = [a for a in range(1, k) if gcd(a, k) == 1]
    
    # Step 3: Calculate N_a = gcd(k, a-1) - 1 for each residue class
    # Note: For p congruent to a (mod k), gcd(k, p-1) = gcd(k, a-1).
    # For a=1, p-1 is a multiple of k, so gcd(k, p-1)=k. We use gcd(k,0)=k.
    sum_N_a = 0
    print("Calculating the exponent w:")
    print(f"k = {k}, phi(k) = {phi_k}")
    print("Residue classes (a mod 12):", residues)
    print("-" * 30)
    
    for a in residues:
        # The value of gcd(12, p-1) for p = a (mod 12) is gcd(12, a-1).
        # math.gcd(12, 0) correctly returns 12.
        g = gcd(k, a - 1)
        N_a = g - 1
        sum_N_a += N_a
        print(f"For a = {a:2d}: N_a = gcd(12, {a-1:2d}) - 1 = {g:2d} - 1 = {N_a}")
        
    print("-" * 30)
    print(f"Sum of N_a values = {sum_N_a}")
    
    # Step 4: Calculate w
    w = sum_N_a / phi_k
    print(f"w = sum(N_a) / phi(12) = {sum_N_a} / {phi_k} = {w}")
    
    # Step 5: Determine alpha and beta
    alpha = 1
    beta = w - 1
    
    print(f"\nThe asymptotic formula is proportional to X^alpha * log(X)^beta.")
    print(f"alpha = {alpha}")
    print(f"beta = w - 1 = {w} - 1 = {int(beta)}")
    
    # Step 6: Calculate the final sum
    result = alpha + beta
    print("\nThe sum alpha + beta is:")
    print(f"{alpha} + {int(beta)} = {int(result)}")
    
    return int(result)

if __name__ == "__main__":
    final_answer = solve_problem()
    # The final answer is returned in the specified format.
    # print(f"\n<<< {final_answer} >>>")