def solve_knuth_bendix_completion():
    """
    This script explains the step-by-step derivation of new rules using the
    Knuth-Bendix completion algorithm for the given system.
    """

    explanation = """
**Step 1: Orient the Initial Rules**

The term ordering is the Lexicographic Path Ordering (LPO) with precedence `f < g < h`. We use this to orient the given pairs into rules.

1.  Pair `(f(g(x), h(x)), g(x))`:
    - `f(g(x), h(x)) > g(x)` because one of its arguments, `h(x)`, is greater than `g(x)` (since `h > g` in the precedence).
    - This gives us **Rule 1: `f(g(x), h(x)) -> g(x)`**

2.  Pair `(f(y, y), g(h(y)))`:
    - `g(h(y)) > f(y, y)` because `g > f` and `g(h(y))` is greater than `y` (which is an argument of `f`).
    - This gives us **Rule 2: `g(h(y)) -> f(y, y)`** (Note the orientation is reversed from the input pair).

3.  Pair `(f(g(x), h(y)), h(x))`:
    - These two terms are incomparable under the given LPO. Neither is consistently greater than the other.
    - Therefore, this pair cannot form a valid rule and is discarded.

Our initial rule set `R` for the completion algorithm is:
- R1: `f(g(x), h(x)) -> g(x)`
- R2: `g(h(y)) -> f(y, y)`

**Step 2: Find Critical Pairs**

We search for overlaps between the left-hand sides (LHS) of the rules in `R`. The only non-trivial overlap is between R2 and R1.
- We unify the LHS of R2, `g(h(y))`, with the subterm `g(x)` in the LHS of R1, `f(g(x), h(x))`.
- The most general unifier is `sigma = {x := h(y)}`.
- The critical pair `(t1, t2)` is formed:
    t1 = sigma applied to `f(r2, h(x))`, which results in `f(f(y,y), h(h(y)))`
    t2 = sigma applied to `r1`, which results in `g(h(y))`
- The critical pair is `(f(f(y,y), h(h(y))), g(h(y)))`.

**Step 3: Resolve the Critical Pair**

We normalize both terms of the pair using the current rule set `R`.
- `g(h(y))` is the LHS of R2, so it normalizes to `f(y,y)`.
- `f(f(y,y), h(h(y)))` cannot be reduced by R1 or R2, so it is already in normal form.
- This gives the divergent critical pair: `(f(f(y,y), h(h(y))), f(y,y))`.

**Step 4: Create a New Rule**

We must orient this new pair using LPO. `f(f(y,y), h(h(y))) > f(y,y)` because they have the same head symbol `f`, and the first argument of the left term `f(y,y)` is greater than the first argument of the right term `y`.
- This creates a new rule which we add to our set:
- **New Rule: `f(f(y, y), h(h(y))) -> f(y, y)`**

**Step 5: Check for More Critical Pairs**

A check for new critical pairs between the new rule and the existing rules reveals no further unifiable overlaps. The algorithm terminates.

**Conclusion**

The only rule added by the Knuth-Bendix completion procedure is the one derived above. Since there is only one new rule, no further ordering is necessary.
"""
    print(explanation)
    final_rule = "f(f(y, y), h(h(y))) -> f(y, y)"
    print("Final added rule:")
    print(final_rule)
    print(f"<<<{final_rule}>>>")

solve_knuth_bendix_completion()