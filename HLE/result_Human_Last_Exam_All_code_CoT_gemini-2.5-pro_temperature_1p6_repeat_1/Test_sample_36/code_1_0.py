import math

def get_divisors(n):
    """
    Helper function to get all divisors of a number.
    """
    divs = set()
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n // i)
    return sorted(list(divs))

print("Step 1: Identify the finite Abelian groups of order 18 and their exponents.")
# The structure theorem for finite Abelian groups gives us the groups of order 18 = 2 * 3^2.
# They are Z_2 x Z_9 (which is isomorphic to Z_18) and Z_2 x Z_3 x Z_3.
# Let's call them G1 and G2.
# G1 = Z_18
# G2 = Z_2 x Z_3 x Z_3
print("The two non-isomorphic Abelian groups of order 18 are G1 = Z_18 and G2 = Z_2 x Z_3 x Z_3.")
exp1 = 18
exp2 = math.lcm(2, 3, 3)
print(f"The exponent of G1 (Z_18) is {exp1}.")
print(f"The exponent of G2 (Z_2 x Z_3 x Z_3) is lcm(2, 3, 3) = {exp2}.")

print("\nStep 2: Determine the possible sets S(rho) for each group.")
print("For a finite Abelian group G, the set of eigenvalues S(rho) is always a group of d-th roots of unity (U_d),")
print("where d is a divisor of the exponent of G. The set of all possible unique sets S(rho) for a given group G")
print("is precisely {U_d | d divides exp(G)}.")

print("\nStep 3: Find the divisors for each exponent.")
divs1 = get_divisors(exp1)
divs2 = get_divisors(exp2)
print(f"The divisors of exp(G1) = {exp1} are: {divs1}")
print(f"The divisors of exp(G2) = {exp2} are: {divs2}")

print("\nStep 4: Find the total number of unique sets by combining the results.")
# The unique sets correspond to the unique divisors found.
# We need to find the size of the union of the two sets of divisors.
union_of_divs = sorted(list(set(divs1) | set(divs2)))
num_divs1 = len(divs1)
num_divs2 = len(divs2)
intersection_divs = sorted(list(set(divs1) & set(divs2)))
num_intersection = len(intersection_divs)
result = len(union_of_divs)

print("The unique sets S(rho) from G1 correspond to the divisors: {" + ", ".join(map(str, divs1)) + "}")
print("The unique sets S(rho) from G2 correspond to the divisors: {" + ", ".join(map(str, divs2)) + "}")
print("\nFinal Calculation using the Principle of Inclusion-Exclusion:")
print(f"Number of sets from G1 = |{{{', '.join(map(str, divs1))}}}| = {num_divs1}")
print(f"Number of sets from G2 = |{{{', '.join(map(str, divs2))}}}| = {num_divs2}")
print(f"Number of common sets = |{{{', '.join(map(str, intersection_divs))}}}| = {num_intersection}")
print(f"Total number of unique sets = (Sets from G1) + (Sets from G2) - (Common Sets)")
print(f"Total number of unique sets = {num_divs1} + {num_divs2} - {num_intersection} = {result}")
