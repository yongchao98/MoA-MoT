import itertools

def is_derangement(perm):
    """
    Checks if a given permutation is a derangement.
    A permutation is a derangement if perm[i] != i+1 for all i.
    Assumes perm is a 0-indexed tuple of 1-indexed values.
    e.g., (2, 3, 1) for sigma(1)=2, sigma(2)=3, sigma(3)=1.
    """
    for i, val in enumerate(perm):
        if i + 1 == val:
            return False
    return True

def count_excedances(perm):
    """
    Counts the number of excedances in a permutation.
    An excedance is an index i such that perm[i-1] > i.
    Assumes perm is a 0-indexed tuple of 1-indexed values.
    """
    count = 0
    for i, val in enumerate(perm):
        if val > i + 1:
            count += 1
    return count

def get_derangements_with_excedances(n):
    """
    Finds all derangements for S_n and their corresponding excedance counts.
    """
    if n < 0:
        return []
    
    # Generate all permutations of (1, 2, ..., n)
    s_n = itertools.permutations(range(1, n + 1))
    
    derangements = []
    for p in s_n:
        if is_derangement(p):
            exc_count = count_excedances(p)
            derangements.append({'perm': p, 'exc': exc_count})
            
    return derangements

def solve_problem():
    """
    Solves the three parts of the problem and prints the final answer.
    """
    # (a) Confirm whether H(U_{n-1, E})(t) = t^{n-1} d_n(t)
    # Let's test for n=4 for a non-trivial example.
    n_a = 4
    # Degree of H(U_{n-1, E})(t) is the rank, which is n-1.
    deg_H = n_a - 1
    # Find derangements and their excedances to compute degree of d_n(t)
    derangements_a = get_derangements_with_excedances(n_a)
    if not derangements_a:
      deg_dn = -1 # Degree of zero polynomial
    else:
      deg_dn = max(d['exc'] for d in derangements_a)
    # Degree of the RHS polynomial P(t) = t^(n-1)d_n(t)
    deg_P = (n_a - 1) + deg_dn
    
    # If the degrees don't match, the identity is false.
    answer_a = "No" if deg_H != deg_P else "Yes"

    # (b) State if the leading coefficient of d_n(t) for any n >= 2 is always 1.
    # The leading coefficient is the number of derangements with maximum excedances.
    is_always_one = True
    for n_b in range(2, 6): # Check for several n
        derangements_b = get_derangements_with_excedances(n_b)
        if not derangements_b:
            is_always_one = False
            break
        max_exc = max(d['exc'] for d in derangements_b)
        leading_coeff = sum(1 for d in derangements_b if d['exc'] == max_exc)
        if leading_coeff != 1:
            is_always_one = False
            break
    answer_b = "Yes" if is_always_one else "No"
    
    # (c) Give the value of d_3(1)
    # d_n(1) is the number of derangements of n elements.
    n_c = 3
    num_derangements_c = len(get_derangements_with_excedances(n_c))
    answer_c = num_derangements_c
    
    # Print the answer in the required format
    print(f"{answer_a}; {answer_b}; {answer_c}")

# Note: The prompt asks me to act as an AI assistant solving a task.
# So I should directly call the function and print the final answer block.
# I'll create the final formatted answer string as the prompt requested.

final_answer_string = ""
# Part a
n_a = 4
deg_H = n_a - 1
derangements_a = get_derangements_with_excedances(n_a)
deg_dn = max(d['exc'] for d in derangements_a)
deg_P = (n_a - 1) + deg_dn
answer_a = "No" if deg_H != deg_P else "Yes"

# Part b
is_always_one = True
for n_b in range(2, 6):
    derangements_b = get_derangements_with_excedances(n_b)
    if not derangements_b:
      is_always_one = False
      break
    max_exc = max(d['exc'] for d in derangements_b)
    leading_coeff = sum(1 for d in derangements_b if d['exc'] == max_exc)
    if leading_coeff != 1:
        is_always_one = False
        break
answer_b = "Yes" if is_always_one else "No"

# Part c
n_c = 3
answer_c = len(get_derangements_with_excedances(n_c))

print(f"<<<{answer_a}; {answer_b}; {answer_c}>>>")