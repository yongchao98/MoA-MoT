import math

def calculate_rigid_rank(N, delta):
    """
    Calculates the largest possible rank r for a (delta, r)-rigid matrix of size NxN
    based on known existence bounds.

    The existence of a (delta, r)-rigid matrix is proven if r satisfies:
    r * (2N - r) / N^2 < 1 - H2(delta)
    where H2(delta) is the binary entropy function.

    This function calculates the upper bound for r based on this inequality.
    """
    print(f"Given parameters: N = {N}, delta = {delta}\n")

    if not (0 < delta < 0.5):
        print("Error: The parameter 'delta' must be in the range (0, 0.5) for the entropy bound to be meaningful.")
        return

    # Step 1: Calculate the binary entropy H2(delta)
    try:
        h2_delta = -delta * math.log2(delta) - (1 - delta) * math.log2(1 - delta)
        print(f"Step 1: Calculate binary entropy H2(delta)")
        print(f"H2({delta}) = -{delta}*log2({delta}) - {1-delta}*log2({1-delta}) = {h2_delta:.4f}")
    except ValueError:
        print("Error: math domain error in entropy calculation. Delta must be between 0 and 1.")
        return
        
    if h2_delta > 1:
        print("Entropy is greater than 1, which implies no rigid matrices can be guaranteed by this method.")
        return

    # Step 2: Calculate the constant c derived from the inequality
    # c is the solution for alpha in alpha(2-alpha) < 1 - H2(delta)
    # The boundary is alpha = 1 - sqrt(H2(delta))
    constant_c = 1 - math.sqrt(h2_delta)
    print(f"\nStep 2: Calculate the factor c for the rank r = c * N")
    print(f"c = 1 - sqrt(H2({delta})) = 1 - sqrt({h2_delta:.4f}) = {constant_c:.4f}")

    # Step 3: Calculate the rank r
    max_r_float = constant_c * N
    # The actual rank must be strictly less than this value, so we take the floor.
    max_r = math.floor(max_r_float)

    print(f"\nStep 3: Calculate the largest rank r")
    print(f"r = c * N = {constant_c:.4f} * {N} = {max_r_float:.2f}")

    print("\n-------------------------------------------------------------")
    print(f"The largest rank r for which a ({delta}, r)-rigid {N}x{N} matrix")
    print("can be constructed by an FNP algorithm is slightly less than the calculated bound.")
    print(f"Thus, the largest integer rank r is {max_r}.")
    print("-------------------------------------------------------------\n")

if __name__ == '__main__':
    # Example values for N and delta
    N_val = 1000
    delta_val = 0.05
    calculate_rigid_rank(N_val, delta_val)