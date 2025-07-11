from z3 import *

def solve_logic_puzzle():
    """
    Solves the logic puzzle by translating sentences into Z3's logic and checking for entailment or contradiction.
    """
    # 1. Setup
    # Using a DeclareSort for Person is more robust for quantifiers
    # than assuming a fixed number of people.
    Person = DeclareSort('Person')
    s = Solver()
    s.set("timeout", 60000)  # 60-second timeout to prevent infinite loops

    # Define predicates for all attributes
    attributes = [
        'Brave', 'Old', 'Curious', 'Creative', 'Calm', 'Quiet', 'Tall',
        'Patient', 'Humble', 'Wise', 'Generous', 'Happy', 'Funny', 'Strong', 'Kind'
    ]
    preds = {name: Function(name, Person, Bool) for name in attributes}

    # Location predicate
    InRoom = Function('InRoom', Person, Bool)
    preds['InRoom'] = InRoom

    # Relations between people
    Older = Function('Older', Person, Person, Bool)
    Richer = Function('Richer', Person, Person, Bool)

    # Declare variables for quantifiers
    x, y = Consts('x y', Person)

    # Unpack predicates into individual variables for cleaner code
    Brave, Old, Curious, Creative, Calm, Quiet, Tall, Patient, Humble, Wise, Generous, Happy, Funny, Strong, Kind = \
        [preds[name] for name in attributes]
    InRoom = preds['InRoom']

    # 2. Assert basic facts for a non-empty universe
    # There must be at least one person in the room and one outside.
    s.add(Exists(x, InRoom(x)))
    s.add(Exists(x, Not(InRoom(x))))

    # 3. Translate and add all 16 premises to the solver
    
    # Premise 1: "A unless B" is equivalent to "B or A"
    p1_A = ForAll(x, Implies(Or(Not(Brave(x)), Old(x)), Not(Curious(x))))
    p1_B = ForAll(x, Implies(Not(InRoom(x)), Not(Creative(x))))
    s.add(Or(p1_B, p1_A))

    # Premise 2: "if C then A else B" is If(C, A, B). Can be simplified.
    p2_C = ForAll(x, Implies(InRoom(x), And(Not(Calm(x)), Not(Brave(x)))))
    p2_B = ForAll(x, Implies(Or(Quiet(x), Not(Creative(x))), Calm(x)))
    s.add(Or(p2_C, p2_B)) # Simplified from If(p2_C, p2_A, p2_B)

    # Premise 3: "A unless B" is "B or A". "iff and vice versa" is "==".
    p3_A = ForAll(x, Implies(InRoom(x), Old(x) == Not(Quiet(x))))
    p3_B = ForAll(x, Implies(InRoom(x), And(Not(Tall(x)), Not(Quiet(x)))))
    s.add(Or(p3_B, p3_A))

    # Premise 4: "if A then B else C" is If(A, B, C)
    p4_A = Not(Exists(x, And(InRoom(x), Curious(x), Wise(x), Not(Tall(x)))))
    p4_B = ForAll(x, Implies(InRoom(x), Implies(Humble(x), Not(Patient(x)))))
    p4_C = ForAll(x, Implies(InRoom(x), Wise(x)))
    s.add(If(p4_A, p4_B, p4_C))

    # Premise 5: "if A then (contradiction) else C" simplifies to "~A and C"
    s.add(Exists(x, Generous(x)))  # not A
    s.add(Exists([x, y], And(Curious(x), Brave(x), Funny(x), Or(Quiet(y), Not(Creative(y))), Older(x, y))))  # C

    # Premise 6: "if A then B else C". "either...or...but not both" is Xor.
    p6_A = ForAll(x, Implies(InRoom(x), Xor(Creative(x), Or(Not(Tall(x)), Not(Generous(x))))))
    p6_B = ForAll(x, Implies(InRoom(x), And(Not(Brave(x)), Creative(x))))
    p6_C = ForAll(x, Implies(InRoom(x), Implies(Patient(x), Not(Wise(x)))))
    s.add(If(p6_A, p6_B, p6_C))

    # Premise 7: "A if B and vice versa" is A == B
    p7_A = ForAll(x, Implies(InRoom(x), And(Not(Patient(x)), Kind(x))))
    p7_B = ForAll(x, Implies(InRoom(x), Generous(x)))
    s.add(p7_A == p7_B)

    # Premise 8: "A unless B" is "B or A"
    p8_A = ForAll(x, And(Generous(x), Not(Quiet(x)), Not(Kind(x))))
    p8_B = ForAll(x, Implies(InRoom(x), Generous(x)))
    s.add(Or(p8_B, p8_A))

    # Premise 9: "A unless B" is "B or A"
    p9_A = Exists([x, y], And(Not(Tall(x)), Not(Strong(x)), Not(Brave(x)), Creative(y), Curious(y), Richer(x, y)))
    p9_B = ForAll(x, Implies(InRoom(x), Xor(Not(Kind(x)), Not(Strong(x)))))
    s.add(Or(p9_A, p9_B))

    # Premise 10: "if A then B else C"
    p10_A = ForAll(x, Implies(InRoom(x), And(Wise(x), Old(x))))
    p10_B = ForAll(x, Implies(InRoom(x), Implies(Calm(x), Or(Not(And(Generous(x), Happy(x))), Not(Wise(x))))))
    p10_C = ForAll(x, Implies(InRoom(x), Not(Generous(x))))
    s.add(If(p10_A, p10_B, p10_C))
    
    # Premise 11
    s.add(ForAll(x, Implies(And(Not(Quiet(x)), Happy(x)), Or(Curious(x), Not(Tall(x))))))

    # Premise 12
    s.add(ForAll(x, Implies(Strong(x), Not(Wise(x)))))

    # Premise 13: "A only if B" is A -> B
    s.add(ForAll(x, Implies(InRoom(x), Implies(And(Wise(x), Not(Humble(x))), And(Not(Quiet(x)), Calm(x))))))

    # Premise 14: The condition contains "Brave and not Brave", which is always false. "False -> A" is always true (a tautology). It provides no new information, so it's ignored.

    # Premise 15
    s.add(ForAll(x, Implies(Not(InRoom(x)), And(Calm(x), Creative(x), Brave(x)))))

    # Premise 16
    s.add(ForAll(x, Implies(InRoom(x), Not(Funny(x)))))

    # 4. Define the Proposition to be tested
    proposition = ForAll(x, Implies(InRoom(x), Tall(x)))

    # 5. Perform the logical checks
    
    # Check for Paradox (self-contradiction in premises)
    if s.check() == unsat:
        print("Analysis Result: The premises are self-contradictory.")
        print("This means the set of statements is a paradox.")
        print("Final Answer: Paradox (uses all premises)")
        return "C"

    # Check for Entailment
    s.push()  # Create a checkpoint
    s.add(Not(proposition))  # Add the negation of the proposition
    entailment_check_result = s.check()
    s.pop()  # Restore to the checkpoint

    # Check for Contradiction
    s.push()
    s.add(proposition) # Add the proposition itself
    contradiction_check_result = s.check()
    s.pop()

    # 6. Interpret the results
    if entailment_check_result == unsat:
        print("Analysis Result: The premises and the negation of the proposition are contradictory.")
        print("This means the premises logically entail the proposition.")
        print("Final Answer: Entailment (uses all premises)")
        return "G"
    elif contradiction_check_result == unsat:
        print("Analysis Result: The premises and the proposition are contradictory.")
        print("This means the premises logically contradict the proposition.")
        print("Final Answer: Contradiction (uses all premises)")
        return "D"
    else:
        print("Analysis Result: The proposition is consistent with the premises, and so is its negation.")
        print("This means the proposition is neutral.")
        print("Final Answer: Neutral")
        return "A"

if __name__ == '__main__':
    final_answer_code = solve_logic_puzzle()
    # The final output format expected by the system.
    # print(f'<<<{final_answer_code}>>>')