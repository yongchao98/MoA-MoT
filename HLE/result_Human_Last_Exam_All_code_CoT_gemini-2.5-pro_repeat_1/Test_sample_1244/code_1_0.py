import math

def gcd_list(numbers):
    """Computes the GCD of a list of integers."""
    if not numbers:
        return 0
    # The GCD of a set of numbers is the same as the GCD of their absolute values.
    numbers = [abs(n) for n in numbers if n != 0]
    if not numbers:
        return 0
    result = numbers[0]
    for i in range(1, len(numbers)):
        result = math.gcd(result, numbers[i])
    return result

def find_vector_for_b():
    """
    This function verifies the claim in part (b) by finding a suitable vector.
    It searches for a primitive vector v in Z^14 such that v.v is a multiple of 54.
    """
    n = 14
    # We search for a vector with a small number of non-zero components
    # to find a solution efficiently. We test for 3 non-zero components.
    search_limit = 10
    for v1 in range(1, search_limit):
        for v2 in range(v1 + 1): # Optimization to avoid permutations
            for v3 in range(v2 + 1):
                components = [v1, v2, v3]
                
                v = components + [0] * (n - len(components))
                
                # Check for primitivity
                if gcd_list(components) != 1:
                    continue
                
                norm_sq = sum(c*c for c in v)
                
                # Check if norm_sq is a multiple of 54
                if norm_sq > 0 and norm_sq % 54 == 0:
                    x_norm_sq = norm_sq / 9
                    print(f"Verification for (b):")
                    print(f"Found a primitive vector v in Z^{n}: {components + [0]*(n-len(components))}")
                    print(f"Primitivity check: gcd{tuple(components)} = {gcd_list(components)}")
                    print(f"The squared norm of v is v.v = {' + '.join([str(c**2) for c in components])} = {norm_sq}")
                    print(f"This norm is a multiple of 54, since {norm_sq} / 54 = {norm_sq // 54}.")
                    print(f"This allows for the construction of an integral 3-neighbor lattice L.")
                    print(f"The corresponding vector x in L has squared norm x.x = (v.v)/9 = {norm_sq} / 9 = {int(x_norm_sq)}")
                    print(f"Checking the final condition: x.x mod 6 = {int(x_norm_sq)} mod 6 = {int(x_norm_sq) % 6}")
                    print(f"Since the result is 0, the condition x.x = 0 (mod 6) is satisfied.")
                    return

# Run the verification function
find_vector_for_b()