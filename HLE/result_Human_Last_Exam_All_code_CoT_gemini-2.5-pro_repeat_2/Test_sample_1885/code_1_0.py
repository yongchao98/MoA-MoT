import sys

# On Python 2, 'print' is a statement, on Python 3, it's a function.
# This is to make the code compatible with both.
try:
    # try to use Python 3's print function
    from __future__ import print_function
except ImportError:
    # Python 2.x
    pass

def explain_problem_setup():
    """Prints an explanation of the problem setup."""
    print("Step 1: Understanding the premises.")
    print("We are given a sequence of functions <f_alpha : alpha < omega_2>.")
    print("Each function f_alpha maps omega_1 to omega_1.")
    print("The sequence is increasing modulo finite, meaning for alpha < beta,")
    print("the set {gamma in omega_1 : f_beta(gamma) <= f_alpha(gamma)} is finite.")
    print("The question is: Must there exist an uncountable set X subset of omega_2 and a function g")
    print("such that for all beta in X, f_beta(gamma) < g(gamma) for all gamma in omega_1?")
    print("-" * 20)

def introduce_coloring_function():
    """Step 2: Define a coloring function based on the sequence."""
    print("Step 2: Defining a coloring function.")
    print("For any pair of ordinals alpha < beta < omega_2, the set where f_beta is not greater than f_alpha is finite.")
    print("Let E_{alpha, beta} = {gamma in omega_1 : f_beta(gamma) <= f_alpha(gamma)}.")
    print("Since E_{alpha, beta} is a finite set of countable ordinals, its supremum is a countable ordinal.")
    print("Let eta_{alpha, beta} = sup(E_{alpha, beta}) + 1. If E is empty, let it be 0.")
    print("So, eta_{alpha, beta} is a countable ordinal (i.e., eta_{alpha, beta} < omega_1).")
    print("For all gamma >= eta_{alpha, beta}, we have f_alpha(gamma) < f_beta(gamma).")
    print("This defines a coloring function c({alpha, beta}) = eta_{alpha, beta} mapping pairs from omega_2 to an ordinal in omega_1.")
    print("-" * 20)

def apply_hajnal_partition_theorem():
    """Step 3: Apply a combinatorial theorem."""
    print("Step 3: Applying Hajnal's Partition Theorem.")
    print("Hajnal's theorem states that for any function c: [omega_2]^2 -> omega_1,")
    print("there exists an uncountable subset H of omega_2 such that the set of colors c([H]^2) is countable.")
    print("Let's apply this to our function eta_{alpha, beta}.")
    print("We get an uncountable set X subset of omega_2 such that the set of values {eta_{alpha, beta} : alpha, beta in X, alpha < beta} is a countable set of ordinals.")
    print("-" * 20)

def analyze_the_monochromatic_set():
    """Step 4: Analyze the properties of the resulting set."""
    print("Step 4: Analyzing the uncountable subfamily.")
    print("Let C = {eta_{alpha, beta} : alpha, beta in X, alpha < beta}. C is a countable set of countable ordinals.")
    print("Let eta_star = sup(C). Since C is countable, eta_star < omega_1.")
    print("This means for any alpha < beta in our set X, eta_{alpha, beta} <= eta_star.")
    print("So, for any alpha < beta in X, and for any gamma > eta_star, we have f_alpha(gamma) < f_beta(gamma).")
    print("-" * 20)

def show_unboundedness():
    """Step 5: Show that this subfamily is unbounded."""
    print("Step 5: Proving the subfamily is unbounded.")
    print("Let's pick an uncountable subset X' of X of size omega_1, and enumerate it as {beta_xi : xi < omega_1}.")
    print("Fix any coordinate gamma > eta_star.")
    print("Consider the sequence of ordinals <f_{beta_xi}(gamma) : xi < omega_1>.")
    print("From the previous step, for any xi < zeta < omega_1, we have f_{beta_xi}(gamma) < f_{beta_zeta}(gamma).")
    print("This is a strictly increasing sequence of ordinals of length omega_1.")
    print("All these ordinals are in omega_1.")
    print("The supremum of such a sequence must be omega_1.")
    print("So, for any gamma > eta_star, sup_{beta in X'} f_beta(gamma) = omega_1.")
    print("-" * 20)

def state_conclusion():
    """Step 6: State the final conclusion."""
    print("Step 6: Conclusion.")
    print("For the family {f_beta : beta in X'} to be bounded by a function g, we would need:")
    print("g(gamma) > sup_{beta in X'} f_beta(gamma) for all gamma.")
    print("But for any gamma > eta_star, this supremum is omega_1.")
    print("So we would need g(gamma) >= omega_1. But g must map to omega_1, so g(gamma) < omega_1.")
    print("This is a contradiction.")
    print("Therefore, the uncountable family {f_beta : beta in X'} is not bounded by any function g.")
    print("This shows that for any such sequence, we can find an uncountable subfamily that is NOT bounded.")
    print("So, it is not necessary that a bounded uncountable subfamily exists.")
    print("-" * 20)
    print("The answer to the question is No.")

def main():
    """Main function to print the argument."""
    explain_problem_setup()
    introduce_coloring_function()
    apply_hajnal_partition_theorem()
    analyze_the_monochromatic_set()
    show_unboundedness()
    state_conclusion()

if __name__ == "__main__":
    main()
