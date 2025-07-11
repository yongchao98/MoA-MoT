import math

def factorial(n):
    """Computes factorial of n."""
    if n == 0:
        return 1
    return math.factorial(n)

def calculate_s():
    """Calculates the total number of distributions S."""
    return factorial(25) // (factorial(5)**5)

def calculate_f():
    """
    Calculates the number of favorable distributions F based on an enumeration
    of symmetric distribution matrices.
    """
    f_fixed = 0
    
    # For each matrix type, we calculate the number of ways W(C) and multiply
    # by the number of such matrices.

    # Type A: diag=(5,5,5,5,5), C=5I. 1 matrix.
    # Hand is (5,0,0,0,0). W_hand = 5!/5! = 1.
    # W(C) = 1^5 = 1.
    w_a = 1 * 1
    f_fixed += w_a

    # Type B: diag=(4,4,5,5,5). C(5,2)=10 matrices.
    # Hands are (4,1,0,0,0) -> 5!/4!=5 and (5,0,0,0,0) -> 1.
    # W(C) = 5 * 5 * 1 * 1 * 1 = 25.
    w_b = 10 * 25
    f_fixed += w_b

    # Type C: diag=(3,3,3,3,3), off-diag is a 5-cycle. (5-1)!/2=12 matrices.
    # Hand is (3,1,1,0,0) -> 5!/(3!1!1!) = 20.
    # W(C) = 20^5.
    w_c = 12 * (20**5)
    f_fixed += w_c

    # Type D: diag=(4,4,4,4,5), off-diag is 2 disjoint edges. 15 matrices.
    # Hands are (4,1,0,0,0) -> 5 and (5,0,0,0,0) -> 1.
    # W(C) = 5^4 * 1 = 625.
    w_d = 15 * 625
    f_fixed += w_d
    
    # Type F: diag=(3,4,4,5,5), off-diag is P3 path. 30 matrices.
    # Hands: (3,1,1,0,0)->20, (4,1,0,0,0)->5, (5,0,0,0,0)->1
    # W(C) = 20^1 * 5^2 * 1^2 = 500.
    w_f = 30 * 500
    f_fixed += w_f

    # Type G: diag=(3,3,4,4,5), off-diag is P4 path. 60 matrices.
    # Hands: (3,1,1,0,0)->20, (4,1,0,0,0)->5, (5,0,0,0,0)->1
    # W(C) = 20^2 * 5^2 * 1^1 = 10000.
    w_g = 60 * 10000
    f_fixed += w_g

    # Type H: diag=(3,3,5,5,5), off-diag has entries of 2. C(5,2)=10 matrices.
    # Hands: (3,2,0,0,0)->5!/(3!2!)=10, (5,0,0,0,0)->1
    # W(C) = 10^2 * 1^3 = 100.
    w_h = 10 * 100
    f_fixed += w_h

    # Type I: diag=(3,3,3,5,5), off-diag is C3 cycle. C(5,3)=10 matrices.
    # Hands: (3,1,1,0,0)->20, (5,0,0,0,0)->1
    # W(C) = 20^3 * 1^2 = 8000.
    w_i = 10 * 8000
    f_fixed += w_i

    # Total F is F_fixed multiplied by 5! for the permutations of dominant types.
    return factorial(5) * f_fixed

def main():
    """
    Main function to calculate and print the probability.
    """
    S = calculate_s()
    F = calculate_f()
    
    # The probability P is F/S.
    # We print the components of the final fraction.
    print(f"The total number of ways to distribute the items is S.")
    print(f"S = 25! / (5!)^5 = {S}")
    print(f"\nThe number of favorable distributions is F.")
    print(f"F = {F}")
    print(f"\nThe probability P = F / S.")
    # To avoid floating point inaccuracies in the final result, we present the fraction.
    # The final result is a fraction, so we simplify it by finding the greatest common divisor.
    common_divisor = math.gcd(F, S)
    f_simplified = F // common_divisor
    s_simplified = S // common_divisor
    
    print(f"P = {F} / {S} = {f_simplified} / {s_simplified}")
    
    # For a decimal approximation:
    # P_decimal = F / S
    # print(f"P ≈ {P_decimal}")

if __name__ == "__main__":
    main()
