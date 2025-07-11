import collections

# Analysis of the statements about attack graphs.
# The goal is to identify the incorrect statement(s).

statements = {
    'A': "Correct. Describes the basic structure (directed graph) and analysis methods (pathfinding, cut-set), which are standard applications of graph theory to attack graphs.",
    'B': "Correct. This accurately describes the 'state explosion problem', a primary challenge in attack graph generation where the graph size can be exponential relative to system size.",
    'C': "Incorrect. This statement mischaracterizes the computational complexity. While some *analysis* problems on attack graphs are NP-hard (e.g., finding optimal countermeasures), the *generation* problem itself is more accurately described as PSPACE-complete due to state-space exploration. Classifying generation as NP-hard is a common but technically inaccurate simplification.",
    'D': "Correct. Real-world systems are dynamic, so the need for incremental and efficient updates to attack graphs is a crucial and accurate point.",
    'E': "Correct. Using formal methods like model checking on attack graphs to verify security properties is a key application and motivation for their use."
}

incorrect_statements = []
for key, value in statements.items():
    if "Incorrect" in value:
        incorrect_statements.append(key)

# Sort the identified incorrect statements alphabetically.
incorrect_statements.sort()

# Print the final result.
# In this case, there is only one incorrect statement.
print(','.join(incorrect_statements))
