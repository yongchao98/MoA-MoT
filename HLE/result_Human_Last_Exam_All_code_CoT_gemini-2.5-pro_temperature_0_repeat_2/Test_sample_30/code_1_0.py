# The user wants to identify the correct statement about interactive proof systems.
# After analyzing all the options, I have determined that option D is the correct one.
# It accurately describes the computational requirements for the prover and the verifier.
# The verifier must be efficient (polynomial-time), while the prover is assumed to be computationally powerful (potentially requiring exponential time or more).

correct_option = 'D'

print(f"The correct statement is D.")
print("Explanation: In an interactive proof system, the verifier is required to be computationally efficient (i.e., run in polynomial time) to be a practical model for verification.")
print("The prover, on the other hand, is assumed to have access to vast computational resources, often modeled as being all-powerful or having at least exponential-time capabilities, in order to find the proof for the given statement.")
print("This asymmetry is a defining characteristic of the model and is what allows it to solve problems in classes like PSPACE, which are believed to be much harder than P.")