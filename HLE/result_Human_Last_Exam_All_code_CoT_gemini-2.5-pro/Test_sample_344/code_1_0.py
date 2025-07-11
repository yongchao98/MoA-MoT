from z3 import *

def solve_logic_puzzle():
    # 1. Define the domain and the predicates (properties)
    Person = DeclareSort('Person')
    
    # Predicates for properties
    properties = [
        "brave", "old", "curious", "creative", "calm", "quiet", "tall",
        "wise", "patient", "humble", "generous", "funny", "happy", "strong", "kind"
    ]
    
    # Functions mapping a Person to a Boolean (True/False)
    p_funcs = {prop: Function(prop.capitalize(), Person, Bool) for prop in properties}
    
    # Predicates for location
    In = Function('In', Person, Bool)
    def Out(p):
        return Not(In(p))

    # Relations between people
    Older = Function('Older', Person, Person, Bool)
    Richer = Function('Richer', Person, Person, Bool)

    # Create a solver instance
    s = Solver()

    # Universal variables for use in quantifiers
    x, y, z = Consts('x y z', Person)

    # Helper to print which premise is being added
    def add_premise(n, p):
        s.add(p)
        print(f"Premise {n} processed.")

    # 2. Translate all premises into Z3 logic
    
    # Premise 1: "if ... then ... unless ..." is ¬B → A or A ∨ B
    p1_A = ForAll(x, Implies(Or(Not(p_funcs["brave"](x)), p_funcs["old"](x)), Not(p_funcs["curious"](x))))
    p1_B = ForAll(x, Implies(Out(x), Not(p_funcs["creative"](x))))
    add_premise(1, Or(p1_A, p1_B))
    
    # Premise 2: "if A then B otherwise C" is (A → B) ∧ (¬A → C)
    p2_A = ForAll(x, Implies(In(x), And(Not(p_funcs["calm"](x)), Not(p_funcs["brave"](x)))))
    p2_B = ForAll(x, Implies(In(x), Not(p_funcs["brave"](x))))
    p2_C = ForAll(x, Implies(Or(p_funcs["quiet"](x), Not(p_funcs["creative"](x))), p_funcs["calm"](x)))
    add_premise(2, And(Implies(p2_A, p2_B), Implies(Not(p2_A), p2_C)))

    # Premise 3: "A unless B" is A ∨ B
    p3_A = ForAll(x, Implies(In(x), p_funcs["old"](x) == Not(p_funcs["quiet"](x))))
    p3_B = ForAll(x, Implies(In(x), And(Not(p_funcs["tall"](x)), Not(p_funcs["quiet"](x)))))
    add_premise(3, Or(p3_A, p3_B))

    # Premise 4: "if A then B otherwise C"
    p4_A = Not(Exists(x, And(In(x), p_funcs["curious"](x), p_funcs["wise"](x), Not(p_funcs["tall"](x)))))
    p4_B = ForAll(x, Implies(In(x), Implies(p_funcs["humble"](x), Not(p_funcs["patient"](x)))))
    p4_C = ForAll(x, Implies(In(x), p_funcs["wise"](x)))
    add_premise(4, And(Implies(p4_A, p4_B), Implies(Not(p4_A), p4_C)))

    # Premise 5: "if A then B else C", where B is a contradiction (is calm and not calm). Simplifies to ¬A ∧ C
    p5_A = ForAll(x, Not(p_funcs["generous"](x)))
    p5_C = Exists([x, y], And(p_funcs["curious"](x), p_funcs["brave"](x), p_funcs["funny"](x), Or(p_funcs["quiet"](y), Not(p_funcs["creative"](y))), Older(x, y)))
    add_premise(5, And(Not(p5_A), p5_C))

    # Premise 6: "if A then B otherwise C"
    p6_A = ForAll(x, Implies(In(x), Xor(p_funcs["creative"](x), And(Not(p_funcs["tall"](x)), Not(p_funcs["generous"](x))))))
    p6_B = ForAll(x, Implies(In(x), And(Not(p_funcs["brave"](x)), p_funcs["creative"](x))))
    p6_C = ForAll(x, Implies(In(x), Implies(p_funcs["patient"](x), Not(p_funcs["wise"](x)))))
    add_premise(6, And(Implies(p6_A, p6_B), Implies(Not(p6_A), p6_C)))

    # Premise 7: "A if B and vice versa" is A ↔ B
    p7_A = ForAll(x, Implies(In(x), And(Not(p_funcs["patient"](x)), p_funcs["kind"](x))))
    p7_B = ForAll(x, Implies(In(x), p_funcs["generous"](x)))
    add_premise(7, p7_A == p7_B)

    # Premise 8: "A unless B" is A ∨ B
    p8_A = ForAll(x, And(p_funcs["generous"](x), Not(p_funcs["quiet"](x)), Not(p_funcs["kind"](x))))
    p8_B = ForAll(x, Implies(In(x), p_funcs["generous"](x)))
    add_premise(8, Or(p8_A, p8_B))

    # Premise 9: "A unless B" is A ∨ B
    p9_A = Exists([x, y], And(Not(p_funcs["tall"](x)), Not(p_funcs["strong"](x)), Not(p_funcs["brave"](x)), p_funcs["creative"](y), p_funcs["curious"](y), Richer(x, y)))
    p9_B = ForAll(x, Implies(In(x), Xor(Not(p_funcs["kind"](x)), Not(p_funcs["strong"](x)))))
    add_premise(9, Or(p9_A, p9_B))
    
    # Premise 10: "if A then B otherwise C"
    p10_A = ForAll(x, Implies(In(x), And(p_funcs["wise"](x), p_funcs["old"](x))))
    p10_B = ForAll(x, Implies(In(x), Implies(p_funcs["calm"](x), Or(Not(And(p_funcs["generous"](x), p_funcs["happy"](x))), Not(p_funcs["wise"](x))))))
    p10_C = ForAll(x, Implies(In(x), Not(p_funcs["generous"](x))))
    add_premise(10, And(Implies(p10_A, p10_B), Implies(Not(p10_A), p10_C)))

    # Premise 11: "if ... then ..."
    add_premise(11, ForAll(x, Implies(And(Not(p_funcs["quiet"](x)), p_funcs["happy"](x)), Or(p_funcs["curious"](x), Not(p_funcs["tall"](x))))))
    
    # Premise 12: "no one ... if ..."
    add_premise(12, ForAll(x, Implies(p_funcs["strong"](x), Not(p_funcs["wise"](x)))))

    # Premise 13: "A only if B" is A → B
    add_premise(13, ForAll(x, Implies(In(x), Implies(And(p_funcs["wise"](x), Not(p_funcs["humble"](x))), And(Not(p_funcs["quiet"](x)), p_funcs["calm"](x))))))

    # Premise 14: Antecedent is a contradiction (brave and not brave), so a tautology. No information.
    print("Premise 14 processed. (Tautology, no assertion added)")

    # Premise 15
    add_premise(15, ForAll(x, Implies(Out(x), And(p_funcs["calm"](x), p_funcs["creative"](x), p_funcs["brave"](x)))))

    # Premise 16
    add_premise(16, ForAll(x, Implies(In(x), Not(p_funcs["funny"](x)))))

    # Make sure the world is not empty to avoid vacuous truths
    s.add(Exists(x, True))

    # 3. Check for Paradox
    print("\nFinal analysis:")
    if s.check() == unsat:
        print("The premises are contradictory, leading to a Paradox.")
        print("Conclusion: Paradox")
        # Find which answer choice corresponds
        # For this problem, C and F are paradoxes. Let's pick C.
        print("<<<C>>>")
        return
    else:
        print("The premises are consistent (not a paradox).")

    # 4. Define the proposition to test
    proposition = ForAll(x, Implies(In(x), p_funcs["tall"](x)))

    # 5. Check for Entailment
    s.push()
    s.add(Not(proposition))
    if s.check() == unsat:
        print("The premises entail the proposition.")
        print("Conclusion: Entailment")
        # G is entailment (uses all premises)
        print("<<<G>>>")
        s.pop()
        return
    s.pop()
    print("The premises do not entail the proposition.")

    # 6. Check for Contradiction
    s.push()
    s.add(proposition)
    if s.check() == unsat:
        print("The premises contradict the proposition.")
        print("Conclusion: Contradiction")
        # D is contradiction (uses all premises)
        print("<<<D>>>")
        s.pop()
        return
    s.pop()
    print("The premises do not contradict the proposition.")

    # 7. If neither, it's Neutral
    print("\nSince the proposition is neither entailed nor contradicted, its status is Neutral.")
    print("<<<A>>>")


if __name__ == '__main__':
    solve_logic_puzzle()
