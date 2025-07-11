import sys

try:
    from z3 import *
except ImportError:
    print("Z3 is not installed. Please install it by running: pip install z3-solver")
    sys.exit()

def solve_logic_problem():
    """
    This function formalizes and solves the given logic puzzle using the Z3 SMT solver.
    """
    # 1. Declare a 'Person' sort and predicates for all properties.
    Person = DeclareSort('Person')
    x, y, z = Consts('x y z', Person)

    InRoom = Function('InRoom', Person, BoolSort())
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
    Funny = Function('Funny', Person, BoolSort())
    Happy = Function('Happy', Person, BoolSort())
    Kind = Function('Kind', Person, BoolSort())
    Strong = Function('Strong', Person, BoolSort())
    Older = Function('Older', Person, Person, BoolSort())
    Richer = Function('Richer', Person, Person, BoolSort())

    # 2. Translate the 16 premises into Z3's logical expressions.
    premises = {}
    
    # P1: “...” unless “...” -> Q | P
    premises['P1'] = Or(
        ForAll(x, Implies(Not(InRoom(x)), Not(Creative(x)))),
        Exists(x, Implies(Or(Not(Brave(x)), Old(x)), Not(Curious(x))))
    )
    # P2: if P then Q else R -> (P -> Q) & (~P -> R)
    p2_P = ForAll(x, Implies(InRoom(x), And(Not(Calm(x)), Not(Brave(x)))))
    p2_Q = ForAll(x, Implies(InRoom(x), Not(Brave(x))))
    p2_R = ForAll(x, Implies(Or(Quiet(x), Not(Creative(x))), Calm(x)))
    premises['P2'] = And(Implies(p2_P, p2_Q), Implies(Not(p2_P), p2_R))
    # P3: “...” unless “...” -> Q | P
    p3_P = ForAll(x, Implies(InRoom(x), Old(x) == Not(Quiet(x))))
    p3_Q = ForAll(x, Implies(InRoom(x), And(Not(Tall(x)), Not(Quiet(x)))))
    premises['P3'] = Or(p3_Q, p3_P)
    # P4: if P then Q else R
    p4_P = ForAll(x, Implies(InRoom(x), Not(And(Curious(x), Wise(x), Not(Tall(x))))))
    p4_Q = ForAll(x, Implies(InRoom(x), Implies(Humble(x), Not(Patient(x)))))
    p4_R = ForAll(x, Implies(InRoom(x), Wise(x)))
    premises['P4'] = And(Implies(p4_P, p4_Q), Implies(Not(p4_P), p4_R))
    # P5: if P then Q else R
    p5_P = ForAll(x, Not(Generous(x)))
    p5_Q = ForAll(x, Implies(Not(InRoom(x)), And(Not(Calm(x)), Calm(x), Not(Creative(x)))))
    p5_R = Exists([x, y], And(Curious(x), Brave(x), Funny(x), Or(Quiet(y), Not(Creative(y))), Older(x, y)))
    premises['P5'] = And(Implies(p5_P, p5_Q), Implies(Not(p5_P), p5_R))
    # P6: if P then Q else R
    p6_P = ForAll(x, Implies(InRoom(x), Xor(Creative(x), Not(And(Tall(x), Generous(x))))))
    p6_Q = ForAll(x, Implies(InRoom(x), And(Not(Brave(x)), Creative(x))))
    p6_R = ForAll(x, Implies(InRoom(x), Implies(Patient(x), Not(Wise(x)))))
    premises['P6'] = And(Implies(p6_P, p6_Q), Implies(Not(p6_P), p6_R))
    # P7: “...” if “...” and vice versa -> P <-> Q
    p7_P = ForAll(x, Implies(InRoom(x), And(Not(Patient(x)), Kind(x))))
    p7_Q = ForAll(x, Implies(InRoom(x), Generous(x)))
    premises['P7'] = (p7_P == p7_Q)
    # P8: “...” unless “...” -> Q | P
    p8_P = ForAll(x, And(Generous(x), Not(Quiet(x)), Not(Kind(x))))
    p8_Q = ForAll(x, Implies(InRoom(x), Generous(x)))
    premises['P8'] = Or(p8_Q, p8_P)
    # P9: “...” unless “...” -> Q | P
    p9_P = Exists([x, y], And(Not(Tall(x)), Not(Strong(x)), Not(Brave(x)), Creative(y), Curious(y), Richer(x, y)))
    p9_Q = ForAll(x, Implies(InRoom(x), Xor(Not(Kind(x)), Not(Strong(x)))))
    premises['P9'] = Or(p9_Q, p9_P)
    # P10: if P then Q else R
    p10_P = ForAll(x, Implies(InRoom(x), And(Wise(x), Old(x))))
    p10_Q = ForAll(x, Implies(InRoom(x), Implies(Calm(x), Or(Not(And(Generous(x), Happy(x))), Not(Wise(x))))))
    p10_R = ForAll(x, Implies(InRoom(x), Not(Generous(x))))
    premises['P10'] = And(Implies(p10_P, p10_Q), Implies(Not(p10_P), p10_R))
    # P11
    premises['P11'] = ForAll(x, Implies(And(Not(Quiet(x)), Happy(x)), Or(Curious(x), Not(Tall(x)))))
    # P12
    premises['P12'] = ForAll(x, Implies(Strong(x), Not(Wise(x))))
    # P13: P only if Q -> P -> Q
    premises['P13'] = ForAll(x, Implies(InRoom(x), Implies(And(Wise(x), Not(Humble(x))), And(Not(Quiet(x)), Calm(x)))))
    # P14: Antecedent contains B & ~B, a contradiction. False -> P is a tautology.
    premises['P14'] = ForAll(x, Implies(And(Not(Strong(x)), Brave(x), Not(Brave(x))), Humble(x)))
    # P15
    premises['P15'] = ForAll(x, Implies(Not(InRoom(x)), And(Calm(x), Creative(x), Brave(x))))
    # P16
    premises['P16'] = ForAll(x, Implies(InRoom(x), Not(Funny(x))))

    # Define the proposition to test.
    proposition = ForAll(x, Implies(InRoom(x), Tall(x)))

    # Add assumptions that people exist inside and outside the room to avoid vacuous truths.
    p_in = Const('p_in', Person)
    p_out = Const('p_out', Person)
    domain_assumptions = [InRoom(p_in), Not(InRoom(p_out))]
    premise_names = set(premises.keys())

    # 3a. Check for Paradox
    s_paradox = Solver()
    s_paradox.set("core.unsat", True)
    for name, p in premises.items(): s_paradox.assert_and_track(p, name)
    for i, p in enumerate(domain_assumptions): s_paradox.assert_and_track(p, f"domain_{i}")
    
    if s_paradox.check() == unsat:
        core = {str(c) for c in s_paradox.unsat_core()}.intersection(premise_names)
        print("Analysis Result: Paradox. The premises are self-contradictory.")
        if len(core) == len(premises): return "C"
        else: return "F"

    # 3b. Check for Contradiction
    s_contradiction = Solver()
    s_contradiction.set("core.unsat", True)
    for name, p in premises.items(): s_contradiction.assert_and_track(p, name)
    for i, p in enumerate(domain_assumptions): s_contradiction.assert_and_track(p, f"domain_{i}")
    s_contradiction.assert_and_track(proposition, "proposition")
    
    if s_contradiction.check() == unsat:
        core = {str(c) for c in s_contradiction.unsat_core()}.intersection(premise_names)
        print("Analysis Result: Contradiction.")
        print("The proposition 'everyone in the room is a tall person' contradicts the premises.")
        print(f"A conflict was found using {len(core)} of the {len(premises)} premises: {sorted(list(core))}")
        if len(core) == len(premises):
            return "D" # Contradiction (uses all premises)
        else:
            return "B" # Contradiction (uses not all premises)

    # 3c. Check for Entailment
    s_entailment = Solver()
    s_entailment.set("core.unsat", True)
    for name, p in premises.items(): s_entailment.assert_and_track(p, name)
    for i, p in enumerate(domain_assumptions): s_entailment.assert_and_track(p, f"domain_{i}")
    s_entailment.assert_and_track(Not(proposition), "neg_proposition")

    if s_entailment.check() == unsat:
        core = {str(c) for c in s_entailment.unsat_core()}.intersection(premise_names)
        print("Analysis Result: Entailment.")
        print("The proposition is entailed by the premises.")
        print(f"A proof was found using {len(core)} of the {len(premises)} premises: {sorted(list(core))}")
        if len(core) == len(premises):
            return "G" # Entailment (uses all premises)
        else:
            return "E" # Entailment (uses not all premises)

    # 3d. If none of the above, it's Neutral.
    print("Analysis Result: Neutral.")
    print("The proposition is neither entailed nor contradicted by the premises.")
    return "A"

# Execute the analysis and print the final answer
final_answer = solve_logic_problem()
print(f"\nFinal Answer Code: {final_answer}")
print(f"<<<{final_answer}>>>")
