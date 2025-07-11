# The user wants to identify the correct statement about interactive proof systems.
# Based on the step-by-step analysis, statement D is the most accurate.

# Statement D: If a prover and verifier are required, there are exponential
# requirements on the computational power of the prover, whereas the verifier
# is required to run in polynomial time.

# This statement correctly captures the core idea of an interactive proof system:
# 1. The Verifier must be efficient (modeled as a probabilistic polynomial-time machine)
#    for the "verification" to be considered efficient.
# 2. The Prover must be powerful enough to find the proof. In the model for IP=PSPACE,
#    the prover must be able to solve PSPACE problems, which may require exponential time.
#    The model often assumes an all-powerful (computationally unbounded) prover.
#    Therefore, stating "exponential requirements" is a correct lower bound on the
#    prover's necessary power for interesting problem classes.

# All other statements contain clear factual errors about definitions or major theorems
# in complexity theory (e.g., GNI, IP=PSPACE, MIP=NEXP, amplification).

final_answer = 'D'
print(f"The correct statement is D.")
print("Analysis: This statement correctly describes the fundamental asymmetry in computational resources within an interactive proof system. The verifier must be computationally limited (polynomial-time) for the verification to be considered efficient, while the prover is assumed to be computationally powerful (at least exponential, and often modeled as all-powerful) to be able to find and provide the necessary proofs. All other options contain significant factual inaccuracies regarding the theory of interactive proofs.")
print(f"<<<{final_answer}>>>")