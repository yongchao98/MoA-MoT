import textwrap

def explain_inconsistency():
    """
    Explains which axiom is inconsistent with the given type theory setup.
    """

    print("### Step 1: Analyzing the Core Problem ###")
    explanation = """
    The system described has a feature called 'structural recursion'. This is a way to ensure that functions always terminate. A recursive call is only allowed if its argument is structurally 'smaller' than the original input.

    The problem introduces a very strange rule for what counts as 'smaller' (a 'subterm'):
    'a lambda (λ x. f) is a subterm of X...'

    This rule effectively means that *any function can be considered a subterm of any piece of data*. This breaks the guarantee of termination, as a recursive function could call itself on an argument (a function) that isn't getting any smaller.
    """
    print(textwrap.dedent(explanation))

    print("### Step 2: The Consequence - Intensional Analysis ###")
    explanation = """
    This broken recursion rule allows us to write functions that analyze the internal definition of other functions. This is called 'intensional analysis'. For example, one could (in theory) write a function H(f) that returns 'True' if f is defined as 'λx.x' and 'False' otherwise. Standard type theories are 'extensional' and forbid this.
    """
    print(textwrap.dedent(explanation))

    print("### Step 3: Finding the Clash with an Axiom ###")
    explanation = """
    Now we look for an axiom that clashes with this new power. The conflict arises between an 'extensional' view of functions and our new 'intensional' power.

    *   B. Functional Extensionality (funext): This axiom says if two functions f and g always produce the same output for the same input, then they ARE the same function (f = g). This is the principle of 'extensionality'. `funext` allows us to prove two functions are equal even if their internal definitions are different.

    *   D. Uniqueness of Identity Proofs (UIP): This axiom says that for any two things x and y, there is at most one proof that they are equal. This is a very strong statement about the nature of equality itself. It gives the equality relation a rigid, 'intensional' character, meaning proofs of equality themselves have a unique structure.

    The conflict is now set up: `funext` lets us equate functions with different structures, while our weird recursion rule lets us write a program that can tell them apart. UIP is the final ingredient that turns this conflict into a full-blown contradiction.
    """
    print(textwrap.dedent(explanation))

    print("### Step 4: The Contradiction ###")
    explanation = """
    Here is a sketch of how the contradiction is formed:

    1.  Using Functional Extensionality, one can prove that two differently-defined functions, `f` and `g`, are equal (i.e., we get a proof `p: f = g`). A classic example involves functions on the empty type.
    2.  Using the weird subterm rule, we can write a function `H` that performs intensional analysis, such that `H(f)` evaluates to `True` and `H(g)` evaluates to `False`.
    3.  Because we have a proof `p: f = g`, we can substitute `g` for `f` and show that `H(f) = H(g)`.
    4.  This leads to `True = False`, a contradiction that makes the whole system inconsistent.

    Uniqueness of Identity Proofs (UIP) is crucial for making this argument formal and watertight within the type theory. It creates a rigid structure for equality proofs that can be broken by the `H` function. The inconsistency of this specific combination is a known result in type theory, first shown by Martin Hofmann.
    """
    print(textwrap.dedent(explanation))

    print("### Conclusion ###")
    answer = "D. Uniqueness of identity proofs"
    print(f"The axiom that is inconsistent with the described system is: {answer}")


explain_inconsistency()
<<<D>>>