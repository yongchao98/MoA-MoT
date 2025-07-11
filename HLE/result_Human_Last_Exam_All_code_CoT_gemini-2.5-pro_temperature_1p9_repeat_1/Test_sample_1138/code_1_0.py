def explain_digital_signature_properties():
    """
    Analyzes the provided statements about digital signature schemes.
    """
    
    print("Analyzing the options based on standard cryptographic principles:")
    print("-" * 60)
    
    print("Initial Analysis of the Premise:")
    print("The question describes ECDSA as an 'existentially forgeable' scheme. This is a contradiction, as ECDSA is designed to be existentially UNFORGEABLE. A 'forgeable' scheme is, by definition, insecure. We will assume the question meant to say 'existentially unforgeable', which is the standard security property.\n")

    # Analysis of Option A
    print("Analysis of A:")
    print("  Statement: For ECDSA, given (m, sig, pk), an adversary can create a new signature sig' for m with no more than negligible probability.")
    print("  Reasoning: This describes strong existential unforgeability. For a secure scheme like ECDSA, it should be computationally infeasible to create a new, valid signature for a message, even if you've already seen one. The probability of doing so is negligible.")
    print("  Conclusion: Statement A is TRUE.\n")

    # Analysis of Option B
    print("Analysis of B:")
    print("  Statement: For ECDSA, given (m, sig, pk), an adversary can recover the secret key sk with no more than negligible probability.")
    print("  Reasoning: The security of the secret key is the most fundamental property of any public-key cryptosystem. If an adversary could recover the secret key, they could forge any signature. The difficulty of this is based on a hard mathematical problem (ECDLP for ECDSA).")
    print("  Conclusion: Statement B is TRUE.\n")
    
    # Analysis of Option C
    print("Analysis of C:")
    print("  Statement: For some schemes, given (m, pk), an adversary can generate a valid signature with NON-negligible probability.")
    print("  Reasoning: This is the definition of a scheme being forgeable (insecure). If we assume the question is about UNFORGEABLE schemes (like ECDSA), then this statement is false by definition.")
    print("  Conclusion: Statement C is FALSE.\n")

    # Analysis of Option D
    print("Analysis of D:")
    print("  Statement: For ALL schemes, given (sig, pk), an adversary can figure out m with no more than negligible probability.")
    print("  Reasoning: This is not universally true. First, some schemes support message recovery. Second, even for schemes without it (like ECDSA), if the message space is small (e.g., 'yes' or 'no'), an attacker can just try verifying the signature for all possible messages. Therefore, this property does not hold for 'all' schemes.")
    print("  Conclusion: Statement D is FALSE.\n")

    # Final Conclusion
    print("-" * 60)
    print("Final Verdict:")
    print("Statements A and B are both true descriptions of a secure, existentially unforgeable signature scheme like ECDSA.")
    print("Therefore, more than one of the options are true.")

explain_digital_signature_properties()
# The final answer is E, as both A and B are correct.
# No equation is present in this problem.
<<<E>>>