import textwrap

def explain_signature_scheme_properties():
    """
    Analyzes the properties of digital signature schemes as described in the user's question.
    """
    # Step 1: Analyze the question's premise.
    # The question presents a contradiction by classifying ECDSA, a secure scheme, as an example of an "existentially forgeable" (insecure) scheme.
    # An "existentially forgeable" scheme is one where an adversary can create a new, valid message-signature pair with non-negligible probability.
    # ECDSA is designed to be "Existentially Unforgeable under a Chosen Message Attack" (EUF-CMA), meaning this is computationally infeasible.
    # We will evaluate each option based on established cryptographic principles.

    print("Analyzing each statement:\n")

    # Step 2: Evaluate statement A.
    analysis_A = """
    A. For ECDSA: Given m, sig, pk, a computationally bounded adversary can create a new, different signature sig' that is verifiable given pk with no more than negligible probability.

    This statement says it is computationally hard to create a second, different signature for a message you already have a signature for. For a secure scheme like ECDSA, this must be true. If an adversary could easily generate new valid signatures for existing messages, the scheme's integrity would be broken. The probability of success is indeed negligible. Therefore, statement A is TRUE.
    """
    print(textwrap.dedent(analysis_A))

    # Step 3: Evaluate statement B.
    analysis_B = """
    B. For ECDSA: Given m, sig, pk, a computationally bounded adversary can recover the secret key sk with no more than negligible probability.

    This statement says that recovering the secret key from public information is computationally hard. The entire security of ECDSA relies on the difficulty of the Elliptic Curve Discrete Logarithm Problem (ECDLP). Recovering the secret key is considered a "total break" of the system and is designed to be infeasible (i.e., possible with only negligible probability). Therefore, statement B is TRUE.
    """
    print(textwrap.dedent(analysis_B))

    # Step 4: Evaluate statement C.
    analysis_C = """
    C. For some existentially forgeable digital signature schemes: Only given m, pk, a computationally bounded adversary can generate sig' that (sig', m) is verifiable against pk with non-negligible probability.

    This statement claims that there exist *some* signature schemes that are broken in a specific way (vulnerable to selective forgery). An "existentially forgeable" scheme is, by definition, an insecure one. It is a fact that insecure or broken signature schemes exist. The statement's claim is weak ("for some"), making it true by the existence of flawed cryptographic algorithms. Therefore, statement C is TRUE.
    """
    print(textwrap.dedent(analysis_C))

    # Step 5: Evaluate statement D.
    analysis_D = """
    D. For all existentially forgeable digital signature schemes: Only given sig, pk, a computationally bounded adversary can figure out m with no more than negligible probability.

    This makes a universal claim ("for all") that message recovery from a signature is hard, even for insecure schemes. We can test this with a counterexample: "textbook" RSA (without hashing), which is known to be existentially forgeable. In textbook RSA, the signature is `sig = m^d mod n`. To recover the message `m` from the signature, anyone with the public key `(e,n)` can compute `m = sig^e mod n`. Message recovery is trivial, not hard. Since a counterexample exists, the universal claim is false. Therefore, statement D is FALSE.
    """
    print(textwrap.dedent(analysis_D))

    # Step 6: Conclude.
    conclusion = """
    Conclusion:
    - Statement A is TRUE.
    - Statement B is TRUE.
    - Statement C is TRUE.
    - Statement D is FALSE.

    Since more than one of the options A-D are true, the correct option is E.
    """
    print(textwrap.dedent(conclusion))

if __name__ == "__main__":
    explain_signature_scheme_properties()