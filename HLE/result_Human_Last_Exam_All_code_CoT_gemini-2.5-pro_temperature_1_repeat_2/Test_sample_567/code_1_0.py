import math

def count_combinations(limit, a):
    """
    Counts the number of non-negative integer pairs (i, j)
    such that i + a*j <= limit.
    This corresponds to the function N_{1,a}(limit).
    """
    count = 0
    j_max = math.floor(limit / a)
    for j in range(j_max + 1):
        i_max = math.floor(limit - a * j)
        count += (i_max + 1)
    return count

def required_combinations(k):
    """
    Calculates the required number of combinations for the ball,
    which is (k+1)(k+2)/2.
    """
    return (k + 1) * (k + 2) // 2

def check_embedding_condition(a):
    """
    Checks if the volume constraint is the only one for a=n^2.
    """
    if a <= 0 or math.sqrt(a) != math.floor(math.sqrt(a)):
        print(f"The value a={a} is not a perfect square.")
        return

    n = int(math.sqrt(a))
    lamb = n  # This is lambda = sqrt(a)

    print(f"Checking for a = {a} and lambda = {lamb}")
    print("The condition is N(k*lambda) >= (k+1)(k+2)/2 for all k>=1.")
    print("For a=n^2, this becomes an equality.")
    print("-" * 50)

    # Check for the first 5 values of k
    for k in range(1, 6):
        limit = k * lamb
        
        # The equation is N(k*lambda) = (k+1)(k+2)/2
        # We calculate both sides.
        N_val = count_combinations(limit, a)
        required_val = required_combinations(k)

        print(f"For k = {k}:")
        print(f"  The equation to check is: N_{{1,{a}}}({k}*{lamb}) = ( {k}+1 )*( {k}+2 )/2")
        print(f"  Calculated left side N_{{1,{a}}}({limit}) = {N_val}")
        print(f"  Calculated right side = {required_val}")
        if N_val == required_val:
            print("  The equality holds.")
        else:
            print("  The equality does NOT hold.")
        print("-" * 50)

# The value of a where the only obstruction becomes the volume constraint.
# This occurs for a = n^2. The first non-trivial case is a = 4.
a_value = 4
check_embedding_condition(a_value)