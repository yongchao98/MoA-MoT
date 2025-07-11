import math

def check_membership(n, c_S):
    """
    Checks if an integer n is in the set S encoded by the real number c_S.
    c_S = sum_{k in S} 3**(-(k+1))
    
    This function mimics the logic that can be built into a first-order formula.
    """
    if not isinstance(n, int) or n < 0:
        print(f"Input n={n} must be a non-negative integer.")
        return False
        
    # These correspond to y1 = 3**(n+1) and y2 = 3**n
    # This step is existentially definable in the language because exponentiation
    # on natural numbers is Diophantine.
    y1 = 3**(n + 1)
    y2 = 3**n
    
    # These correspond to z1 = floor(y1 * c_S) and z2 = floor(y2 * c_S)
    # The floor function is existentially definable.
    z1 = math.floor(y1 * c_S)
    z2 = math.floor(y2 * c_S)
    
    # This is the final check, a simple polynomial equation.
    # The (n+1)-th trit is 1.
    digit = z1 - 3 * z2
    
    print(f"For n = {n}:")
    print(f"  y1 = 3**(n+1) = {y1}")
    print(f"  y2 = 3**n = {y2}")
    print(f"  c_S = {c_S}")
    print(f"  z1 = floor(y1 * c_S) = floor({y1 * c_S:.4f}) = {z1}")
    print(f"  z2 = floor(y2 * c_S) = floor({y2 * c_S:.4f}) = {z2}")
    print(f"  Digit = z1 - 3 * z2 = {z1} - 3 * {z2} = {digit}")
    
    return digit == 1

# --- Example ---
# Let's define an arbitrary set S, e.g., the set of prime numbers less than 20
S = {2, 3, 5, 7, 11, 13, 17, 19}

# Encode S into the parameter c_S
# This real number c_S is the parameter in our formula
c_S = sum(3**(-(n+1)) for n in S)

# Now, let's test membership for a few numbers
print(f"Testing membership for the set S = {S}")
print("-" * 20)
is_in_set = check_membership(5, c_S)
print(f"Is 5 in S? {is_in_set}")
print("-" * 20)
is_in_set = check_membership(6, c_S)
print(f"Is 6 in S? {is_in_set}")
print("-" * 20)
