import sys

# This script explains and derives the formula for the probability of an edge
# in a jointly exchangeable random graph.

print("Step 1: Understanding Jointly Exchangeable Graphs")
print("A random graph is 'jointly exchangeable' if its probability distribution is the same regardless of how we label the nodes. This is a symmetry property.")
print("The Aldous-Hoover theorem provides a way to represent such graphs. A common and intuitive representation is the W-graph (or graphon) model.")
print("-" * 50)

print("Step 2: The W-graph Model")
print("In this model, the graph is generated as follows:")
print("1. For each node i (from 1 to N), we draw a latent variable u_i, independently from a distribution on [0,1].")
print("   Given the prompt, we assume this distribution is Uniform: u_i ~ U(0,1).")
print("2. The existence of a link y_ij between nodes i and j depends on their latent variables via a symmetric function W:[0,1]x[0,1] -> [0,1], called the graphon.")
print("   The conditional probability of a link is: P(y_ij = 1 | u_i, u_j) = W(u_i, u_j)")
print("-" * 50)

print("Step 3: Calculating the Marginal Probability")
print("We want to find the unconditional (or marginal) probability of a link, P(y_ij = 1).")
print("To do this, we average the conditional probability over all possible values of the latent variables u_i and u_j.")
print("Using the law of total probability, this is written as an expectation:")
print("P(y_ij = 1) = E[ P(y_ij = 1 | u_i, u_j) ] = E[ W(u_i, u_j) ]")
print("-" * 50)

print("Step 4: The Final Formula")
print("Since u_i and u_j are independent draws from U(0,1), their joint probability density function is p(u_i, u_j) = 1 on the unit square [0,1]x[0,1].")
print("The expectation is therefore a double integral over this unit square.")
print("\nThe final formula for the probability of a link y_ij is:")
# For a prettier output in a terminal that supports unicode
if sys.stdout.encoding.lower().startswith('utf-8'):
    print(u"P(y_ij = 1) = \u222B\u222B W(u,v) du dv  (integrated over u in [0,1] and v in [0,1])")
else:
    print("P(y_ij = 1) = integral from 0 to 1 (for u) integral from 0 to 1 (for v) of W(u,v) du dv")

print("\nThis value represents the expected density of edges in the graph.")
