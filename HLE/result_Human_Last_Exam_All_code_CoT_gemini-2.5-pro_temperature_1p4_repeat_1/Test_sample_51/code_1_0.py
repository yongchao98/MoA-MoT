import sys

# In dependent type theory, proofs are computational objects. We can model them as classes.
# This class represents a proof of equality, carrying a "reason" for how it was constructed.
class EqualityProof:
    """A toy model for an equality proof object."""
    def __init__(self, term1, term2, reason):
        if term1 != term2:
            # A proof object can only exist if the statement is true.
            raise TypeError("Cannot construct a proof for a false statement.")
        self.term1 = term1
        self.term2 = term2
        self.reason = reason

    def __str__(self):
        return f"Proof({self.term1} == {self.term2}, reason='{self.reason}')"

# --- The Setup ---

# 1. We have canonical proofs from reflexivity. In type theory, this is a constructor.
# Let's create a proof for `1 == 1`.
proof_by_refl = EqualityProof(1, 1, "reflexivity")
print(f"Constructed a proof object: {proof_by_refl}")


# 2. In complex theories, you can construct proofs in other ways.
# Let's imagine a complex theorem `T` also proves `1 == 1`, yielding a different proof object.
# (In a real system, this could be the result of a long proof chain).
proof_by_theorem = EqualityProof(1, 1, "via complex theorem T")
print(f"Constructed another proof object: {proof_by_theorem}")
print("-" * 20)


# 3. "Structural Recursion" on equality proofs is like a function that can inspect the
# structure (the "reason") of the proof object.
def analyze_proof(p: EqualityProof):
    """
    This function mimics pattern matching on the proof's constructor.
    It returns a different value based on how the proof was made.
    """
    print(f"Analyzing proof: {p}...")
    if p.reason == "reflexivity":
        # This branch is taken for proofs constructed "by reflexivity".
        return 1
    else:
        # This branch is taken for all other proofs.
        return 0

# Let's run the analysis on our two proofs.
result1 = analyze_proof(proof_by_refl)
result2 = analyze_proof(proof_by_theorem)

print(f"\nAnalysis of first proof returned: {result1}")
print(f"Analysis of second proof returned: {result2}")
print("-" * 20)

# --- The Inconsistency ---

# 4. Now, let's introduce the axiom of "Uniqueness of Identity Proofs" (UIP).
print("Axiom Added: Uniqueness of Identity Proofs (UIP).")
print("UIP states that any two proofs of the same fact are equal.")
print(f"Under UIP, since both proofs show '1 == 1', we must have:")
print(f"'{proof_by_refl}' == '{proof_by_theorem}'")
print("-" * 20)

# 5. If the proofs are equal, any function applied to them must produce the same result.
print("The consequence of UIP is that `analyze_proof(proof_by_refl)` MUST EQUAL `analyze_proof(proof_by_theorem)`.")

# 6. Here is the contradiction.
print("\n--- The Contradiction ---")
print("Our analysis function showed:")
print(f"analyze_proof(proof_by_refl)  => {result1}")
print(f"analyze_proof(proof_by_theorem) => {result2}")
print("\nBut UIP implies these results must be equal.")
print("Therefore, the system has proven the following equation:")

# Final explicit output of the contradictory equation.
final_lhs = result1
final_rhs = result2
print(f"\n{final_lhs} = {final_rhs}\n")

if final_lhs != final_rhs:
    print("This is a contradiction. The system is inconsistent.")
    print("The inconsistency arises from combining structural recursion (the 'analyze_proof' function) with the UIP axiom.")
else:
    # This case should not be reached with the current logic.
    print("No contradiction found.")
