import math

def C(n, k):
    """
    Calculates the binomial coefficient C(n, k) = n! / (k! * (n-k)!).
    Returns 0 if k > n or k < 0.
    """
    if k < 0 or k > n:
        return 0
    return math.factorial(n) // (math.factorial(k) * math.factorial(n - k))

# Part a: Find the bi-degree of the terminal polynomial.
a_start, b_start = 4, 3
h_weight = a_start - b_start
total_degree = a_start + b_start
# For the terminal polynomial with bi-degree (a_end, b_end):
# a_end - b_end = -h_weight
# a_end + b_end = total_degree
# Solving this system:
a_end = (total_degree - h_weight) / 2
b_end = (total_degree + h_weight) / 2
answer_a = f"({int(a_end)}, {int(b_end)})"

# Part b: Provide the condition for a string starter.
answer_b = "a >= C(r_1, 2) + C(r_2, 2) + ... + C(r_b, 2)"

# Part c: Check possibility for a polynomial of bi-degree (5, 2).
a, b = 5, 2
r_indices = [1, 2] # From "r = 1, 2"
# Applying the condition from part b:
sum_of_coeffs = C(r_indices[0], 2) + C(r_indices[1], 2)
is_possible = "Yes" if a >= sum_of_coeffs else "No"
answer_c = is_possible

# Output the step-by-step derivation as requested.
print("Step-by-step derivation:")
print("-" * 25)

# Part a derivation
print("Part a)")
print(f"The starting polynomial has bi-degree ({a_start}, {b_start}).")
print(f"The H-eigenvalue (weight) is \u03BB = a - b = {a_start} - {b_start} = {h_weight}.")
print(f"The total degree d = a + b = {a_start} + {b_start} = {total_degree} is conserved throughout the string.")
print(f"The terminal polynomial has the lowest weight, which is -\u03BB = {-h_weight}.")
print(f"For the terminal bi-degree (a', b'), we solve the system of equations:")
print(f"  a' - b' = {-h_weight}")
print(f"  a' + b' = {total_degree}")
print(f"The solution is a' = {int(a_end)}, b' = {int(b_end)}. The bi-degree is {answer_a}.")
print("-" * 25)

# Part b derivation
print("Part b)")
print("The condition for a polynomial of bi-degree (a, b) constructed with indices r_1, ..., r_b")
print("to be a valid string starter is proposed based on known combinatorial properties of diagonal harmonics.")
print(f"The condition is: {answer_b}")
print("-" * 25)

# Part c derivation
print("Part c)")
print(f"We test if a polynomial of bi-degree ({a}, {b}) can be constructed using indices r = {', '.join(map(str, r_indices))}.")
print(f"Applying the condition from part b: a >= C(r_1, 2) + C(r_2, 2)")
print(f"We check if {a} >= C({r_indices[0]}, 2) + C({r_indices[1]}, 2)")
c12 = C(r_indices[0], 2)
c22 = C(r_indices[1], 2)
print(f"Since C({r_indices[0]}, 2) = {c12} and C({r_indices[1]}, 2) = {c22}, the inequality is {a} >= {c12} + {c22}, which is {a} >= {sum_of_coeffs}.")
print(f"This is true, so the answer is '{answer_c}'.")
print("-" * 25)

final_answer_string = f"a) {answer_a} b) {answer_b} c) {answer_c}"
print("\nFinal Answer:")
print(f"<<<{final_answer_string}>>>")