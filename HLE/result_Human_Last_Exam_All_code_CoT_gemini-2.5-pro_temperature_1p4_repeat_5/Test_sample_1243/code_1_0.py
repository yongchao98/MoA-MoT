import math

def combinations(n, k):
    """Calculates the binomial coefficient C(n, k)"""
    if k < 0 or k > n:
        return 0
    return math.factorial(n) // (math.factorial(k) * math.factorial(n - k))

print("Calculation Steps:")
print("="*20)

# Part (a)
print("\nPart (a):")
print("We count the number of non-zero binary vectors of length 16 with weight divisible by 4.")
print("The count is C(16, 4) + C(16, 8) + C(16, 12) + C(16, 16).")
c16_4 = combinations(16, 4)
c16_8 = combinations(16, 8)
c16_12 = combinations(16, 12)
c16_16 = combinations(16, 16)
answer_a = c16_4 + c16_8 + c16_12 + c16_16
print(f"The equation is: {c16_4} + {c16_8} + {c16_12} + {c16_16} = {answer_a}")

# Part (b)
print("\nPart (b):")
print("This is a theoretical question. Our analysis shows that the conditions for the neighbor lattice to be even lead to a contradiction with the 3-isotropic property.")
answer_b = "no"
print(f"The answer is {answer_b}.")

# Part (c)
print("\nPart (c):")
print("We count vectors x in Z^12 with norm x.x = 4 that are 2-primitive.")
print("These are the vectors with exactly four components equal to +-1 and others zero.")
print("Number of such vectors = C(12, 4) * 2^4.")
c12_4 = combinations(12, 4)
num_vectors_c = c12_4 * 16
print(f"First, we calculate C(12, 4) = {c12_4}.")
print(f"Number of vectors = {c12_4} * 16 = {num_vectors_c}.")
print("Each neighbor lattice N(x) corresponds to two vectors {x, -x}, so we divide by 2.")
answer_c = num_vectors_c // 2
print(f"Number of lattices = {num_vectors_c} / 2 = {answer_c}.")

print("\n" + "="*20)
print("Final Answer:")
print(f"(a) {answer_a}; (b) {answer_b}; (c) {answer_c}")
