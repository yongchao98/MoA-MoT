def check_unlearnable_mixing_operators():
    """
    Checks which mixing-dimension logical operators are not learnable by an
    additively separable model.

    The model's score function has the form Score = f(p1, p2) + g(q1, q2) + b.
    To learn a target function T(p1, q2), T must also be additively separable.
    The test for additive separability for a binary function T(x, y) is:
    T(1,1) - T(1,0) - T(0,1) + T(0,0) = 0.

    If this interaction term is non-zero, the operator is not learnable.
    """
    
    # Define the mixing-dimension operators as functions
    # T(p1, q2)
    operators = {
        "X' (XOR mixing)": lambda p1, q2: p1 ^ q2,
        "C' (Conjunction mixing)": lambda p1, q2: p1 & q2,
        "D' (Disjunction mixing)": lambda p1, q2: p1 | q2,
        "E' (Equivalence mixing)": lambda p1, q2: 1 if p1 == q2 else 0,
        "I' (Implication mixing)": lambda p1, q2: int((not p1) or q2)
    }

    unlearnable_ops = []
    
    print("Checking which mixing-dimension operators are unlearnable:")
    print("-" * 60)

    for name, op_func in operators.items():
        t11 = op_func(1, 1)
        t10 = op_func(1, 0)
        t01 = op_func(0, 1)
        t00 = op_func(0, 0)
        
        interaction_term = t11 - t10 - t01 + t00
        
        print(f"Operator: {name}")
        # Final Answer requirement: output each number in the final equation!
        print(f"Interaction Term = T(1,1) - T(1,0) - T(0,1) + T(0,0)")
        print(f"                 = {t11} - {t10} - {t01} + {t00} = {interaction_term}")
        
        if interaction_term != 0:
            print("Result: NOT Learnable (interaction term is non-zero)\n")
            unlearnable_ops.append(name.split(" ")[0])
        else:
            # This case won't be hit for the given operators
            print("Result: Potentially Learnable (interaction term is zero)\n")

    print("-" * 60)
    print("Conclusion: All element-wise operators (X, C, D, E, I) are learnable.")
    print("The operators that CANNOT be learned are:")
    print(', '.join(sorted(unlearnable_ops)))

if __name__ == '__main__':
    check_unlearnable_mixing_operators()
