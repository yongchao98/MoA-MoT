# The user is expected to have z3-solver installed (pip install z3-solver)
from z3 import *

def solve_logic_problem():
    """
    This function models and solves the given symbolic logic problem using Z3.
    """
    # 1. Setup
    # Declare a new uninterpreted sort for Person
    Person = DeclareSort('Person')

    # Define the solver
    s = Solver()

    # Define all predicates as functions from Person to Bool
    Brave = Function('Brave', Person, BoolSort())
    Old = Function('Old', Person, BoolSort())
    Curious = Function('Curious', Person, BoolSort())
    Creative = Function('Creative', Person, BoolSort())
    Calm = Function('Calm', Person, BoolSort())
    Quiet = Function('Quiet', Person, BoolSort())
    Tall = Function('Tall', Person, BoolSort())
    Wise = Function('Wise', Person, BoolSort())
    Patient = Function('Patient', Person, BoolSort())
    Humble = Function('Humble', Person, BoolSort())
    Generous = Function('Generous', Person, BoolSort())
    Happy = Function('Happy', Person, BoolSort())
    Funny = Function('Funny', Person, BoolSort())
    Strong = Function('Strong', Person, BoolSort())
    Kind = Function('Kind', Person, BoolSort())

    # Define relations (though they only appear in existentials)
    Older = Function('Older', Person, Person, BoolSort())
    Richer = Function('Richer', Person, Person, BoolSort())

    # Define domains (locations)
    InRoom = Function('InRoom', Person, BoolSort())
    OutsideRoom = Function('OutsideRoom', Person, BoolSort())

    # Quantifier variable
    x, y = Consts('x y', Person)

    # Add axioms about the domains
    # A person is either in the room or outside, but not both.
    s.add(ForAll(x, Xor(InRoom(x), OutsideRoom(x))))
    # Assert that there's at least one person in the room and at least one outside
    # This prevents vacuous truths where a domain is empty.
    p_in = Const('p_in', Person)
    p_out = Const('p_out', Person)
    s.add(InRoom(p_in))
    s.add(OutsideRoom(p_out))

    # 2. Translate Premises into Z3 constraints
    # Helper for 'A unless B' -> Or(B, A)
    # Helper for 'if P then Q otherwise R' -> And(Implies(P,Q), Implies(Not(P),R))
    # Helper for 'A if B and vice versa' -> A == B

    # Premise 1
    p1_A = ForAll(x, Implies(Or(Not(Brave(x)), Old(x)), Not(Curious(x))))
    p1_B = ForAll(x, Implies(OutsideRoom(x), Not(Creative(x))))
    s.add(Or(p1_B, p1_A))

    # Premise 2
    p2_P = ForAll(x, Implies(InRoom(x), And(Not(Calm(x)), Not(Brave(x)))))
    p2_Q = ForAll(x, Implies(InRoom(x), Not(Brave(x)))) # Note: P->Q is a tautology here.
    p2_R = ForAll(x, Implies(Or(Quiet(x), Not(Creative(x))), Calm(x)))
    s.add(And(Implies(p2_P, p2_Q), Implies(Not(p2_P), p2_R)))

    # Premise 3
    p3_A = ForAll(x, Implies(InRoom(x), Old(x) == Not(Quiet(x))))
    p3_B = ForAll(x, Implies(InRoom(x), And(Not(Tall(x)), Not(Quiet(x)))))
    s.add(Or(p3_B, p3_A))

    # Premise 4
    p4_P = Not(Exists(x, And(InRoom(x), Curious(x), Wise(x), Not(Tall(x)))))
    p4_Q = ForAll(x, Implies(InRoom(x), Implies(Humble(x), Not(Patient(x)))))
    p4_R = ForAll(x, Implies(InRoom(x), Wise(x)))
    s.add(And(Implies(p4_P, p4_Q), Implies(Not(p4_P), p4_R)))

    # Premise 5
    p5_P = ForAll(x, Not(Generous(x)))
    p5_Q = ForAll(x, Implies(OutsideRoom(x), And(Not(Calm(x)), Calm(x), Not(Creative(x))))) # Consequent is False
    p5_R = Exists([x, y], And(Curious(x), Brave(x), Funny(x), Or(Quiet(y), Not(Creative(y))), Older(x,y)))
    s.add(And(Implies(p5_P, p5_Q), Implies(Not(p5_P), p5_R)))

    # Premise 6
    p6_P = ForAll(x, Implies(InRoom(x), Xor(Creative(x), And(Not(Tall(x)), Not(Generous(x))))))
    p6_Q = ForAll(x, Implies(InRoom(x), And(Not(Brave(x)), Creative(x))))
    p6_R = ForAll(x, Implies(InRoom(x), Implies(Patient(x), Not(Wise(x)))))
    s.add(And(Implies(p6_P, p6_Q), Implies(Not(p6_P), p6_R)))

    # Premise 7
    p7_A = ForAll(x, Implies(InRoom(x), And(Not(Patient(x)), Kind(x))))
    p7_B = ForAll(x, Implies(InRoom(x), Generous(x)))
    s.add(p7_A == p7_B)

    # Premise 8
    p8_A = ForAll(x, And(Generous(x), Not(Quiet(x)), Not(Kind(x))))
    p8_B = ForAll(x, Implies(InRoom(x), Generous(x)))
    s.add(Or(p8_B, p8_A))

    # Premise 9
    p9_A = Exists([x, y], And(Not(Tall(x)), Not(Strong(x)), Not(Brave(x)), Creative(y), Curious(y), Richer(x,y)))
    p9_B = ForAll(x, Implies(InRoom(x), Xor(Not(Kind(x)), Not(Strong(x)))))
    s.add(Or(p9_B, p9_A))

    # Premise 10
    p10_P = ForAll(x, Implies(InRoom(x), And(Wise(x), Old(x))))
    p10_Q = ForAll(x, Implies(InRoom(x), Implies(Calm(x), Or(Not(And(Generous(x), Happy(x))), Not(Wise(x))))))
    p10_R = ForAll(x, Implies(InRoom(x), Not(Generous(x))))
    s.add(And(Implies(p10_P, p10_Q), Implies(Not(p10_P), p10_R)))

    # Premise 11
    s.add(ForAll(x, Implies(And(Not(Quiet(x)), Happy(x)), Or(Curious(x), Not(Tall(x))))))

    # Premise 12
    s.add(ForAll(x, Implies(Strong(x), Not(Wise(x)))))

    # Premise 13
    s.add(ForAll(x, Implies(InRoom(x), Implies(And(Wise(x), Not(Humble(x))), And(Not(Quiet(x)), Calm(x))))))

    # Premise 14: Antecedent is a contradiction (Brave AND NOT Brave), so statement is always true (tautology)
    # No need to add it, but included for completeness.
    s.add(ForAll(x, Implies(And(Not(Strong(x)), Brave(x), Not(Brave(x))), Humble(x))))
    
    # Premise 15
    s.add(ForAll(x, Implies(OutsideRoom(x), And(Calm(x), Creative(x), Brave(x)))))

    # Premise 16
    s.add(ForAll(x, Implies(InRoom(x), Not(Funny(x)))))

    # 3. Define the Proposition
    proposition = ForAll(x, Implies(InRoom(x), Tall(x)))

    # 4. Analyze
    # Check 1: Are the premises themselves contradictory?
    if s.check() == unsat:
        print("Result: Paradox (uses all premises)")
        print("The set of 16 premises is self-contradictory.")
        return 'C'

    # Check 2: Is the proposition an entailment?
    # (Premises => Prop) is equivalent to (Premises AND NOT Prop) being unsatisfiable
    s.push()
    s.add(Not(proposition))
    is_entailment = (s.check() == unsat)
    s.pop()

    if is_entailment:
        print("Result: Entailment (uses all premises)")
        print("The premises logically entail the proposition. The proposition MUST be true.")
        return 'G'

    # Check 3: Is the proposition a contradiction?
    # (Premises => NOT Prop) is equivalent to (Premises AND Prop) being unsatisfiable
    s.push()
    s.add(proposition)
    is_contradiction = (s.check() == unsat)
    s.pop()

    if is_contradiction:
        print("Result: Contradiction (uses all premises)")
        print("The premises logically contradict the proposition. The proposition MUST be false.")
        return 'D'

    # Check 4: If neither, it's neutral
    print("Result: Neutral")
    print("The proposition is consistent with the premises, but not entailed by them.")
    print("This means worlds can exist where the proposition is true, and worlds can exist where it is false.")
    return 'A'

# Execute the analysis and print the result.
final_answer = solve_logic_problem()
# The strange instruction "output each number in the final equation" doesn't map to a standard logical
# concept. I will interpret it as a request to show which numbered premises were involved in the final
# state. As all premises were used to constrain the model, all numbers from 1 to 16 are relevant.
# To satisfy the instruction symbolically, I'll print the set of numbers.
print("\nRelevant premises for the conclusion: {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16}")

#<<<A>>>
if __name__ == '__main__':
    solve_logic_problem()