import math

def combinations(n, k):
    """Calculates the number of combinations C(n, k)"""
    if k < 0 or k > n:
        return 0
    # Use integer division
    return math.factorial(n) // (math.factorial(k) * math.factorial(n - k))

# Part (a) Calculation
n_a = 16
c_16_4 = combinations(n_a, 4)
c_16_8 = combinations(n_a, 8)
c_16_12 = combinations(n_a, 12)
c_16_16 = combinations(n_a, 16)
answer_a = c_16_4 + c_16_8 + c_16_12 + c_16_16
eq_a = f"C(16,4) + C(16,8) + C(16,12) + C(16,16) = {c_16_4} + {c_16_8} + {c_16_12} + {c_16_16} = {answer_a}"

# Part (b) Answer
answer_b = "no"

# Part (c) Calculation
n_c = 12
answer_c = combinations(n_c, 4)
eq_c = f"C(12,4) = {answer_c}"

# Final Output
print("Calculations:")
print(f"(a) {eq_a}")
print(f"(b) Based on the analysis, it is not possible for the neighbor to be even.")
print(f"(c) {eq_c}")

print("\nFinal Answer:")
print(f"(a) {answer_a}; (b) {answer_b}; (c) {answer_c}")

print("<<<" + f"(a) {answer_a}; (b) {answer_b}; (c) {answer_c}" + ">>>")