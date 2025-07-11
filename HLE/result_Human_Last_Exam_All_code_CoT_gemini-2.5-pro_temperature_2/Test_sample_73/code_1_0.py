import math

def solve_pentagon_genus():
    """
    Calculates the genus of the configuration space of a hinged regular pentagon
    with two adjacent vertices nailed to the plane.
    """
    # The problem concerns a pentagon, so we work with n=5.
    n = 5
    print(f"The problem is for a pentagon, so we set n = {n}.")

    # The genus is found via the Euler characteristic, chi.
    # For a regular n-gon (with n odd), the formula for the Euler characteristic
    # of its configuration space is: chi = (2 - n) * C(n - 3, (n - 3) // 2)
    print("We will use the formula for the Euler characteristic: chi = (2 - n) * C(n - 3, (n - 3) // 2)")

    # Step 1: Calculate the parameters for the binomial coefficient C(k, r)
    k = n - 3
    r = (n - 3) // 2
    print(f"\nThe parameters for the binomial coefficient are k = n - 3 = {k} and r = (n - 3) // 2 = {r}.")

    # Step 2: Calculate the binomial coefficient C(k, r)
    # C(k, r) = k! / (r! * (k - r)!)
    try:
        # math.comb is available in Python 3.8+
        C_k_r = math.comb(k, r)
    except AttributeError:
        # Manual implementation for older Python versions
        C_k_r = math.factorial(k) // (math.factorial(r) * math.factorial(k - r))
    print(f"The binomial coefficient is C({k}, {r}) = {C_k_r}.")

    # Step 3: Calculate the Euler characteristic 'chi'
    chi = (2 - n) * C_k_r
    print("\nCalculating the Euler characteristic 'chi':")
    print(f"chi = ({2} - {n}) * {C_k_r} = {chi}")

    # Step 4: Calculate the genus 'g' from the Euler characteristic
    # The formula relating genus and chi for a smooth orientable surface is: chi = 2 - 2g
    # This can be rearranged to g = (2 - chi) / 2
    g = (2 - chi) // 2
    print("\nThe genus 'g' is related to the Euler characteristic by the formula: chi = 2 - 2g")
    print("Solving for g:")
    print(f"{chi} = 2 - 2g")
    two_g = 2 - chi
    print(f"2g = 2 - ({chi}) = {two_g}")
    print(f"g = {two_g} / 2 = {g}")
    
    print(f"\nThus, the genus of the surface is {g}.")
    
    print("\nVerification of the final equation:")
    print(f"chi = 2 - 2g  =>  {chi} = {2} - {2} * {g}")


solve_pentagon_genus()