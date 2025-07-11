def solve_group_theory_problem():
    """
    This function explains and states the characterization of finite filled nilpotent groups.
    """

    title = "Characterization of Finite Filled Nilpotent Groups"
    
    # Step-by-step reasoning
    reasoning = [
        ("1. Definitions:",
         "A product-free set S in a group G is a subset where for any x, y in S, the product xy is not in S. "
         "A group G is 'filled' if the union of all its maximal product-free sets is equal to G itself. "
         "A finite group is 'nilpotent' if it is the direct product of its Sylow p-subgroups."),
        
        ("2. The Frattini Reduction:",
         "A key theorem by V. A. Petrov states that a finite group G is filled if and only if its Frattini quotient, G/Φ(G), is filled. "
         "(The Frattini subgroup Φ(G) is the intersection of all maximal subgroups of G)."),

        ("3. Nilpotent Group Property:",
         "For any finite nilpotent group G, its Frattini quotient G/Φ(G) is a finite abelian group."),
        
        ("4. Filled Abelian Groups:",
         "Another theorem states that a finite abelian group A is filled if and only if A has a subgroup of index 2 or a subgroup of index 3."),

        ("5. Subgroups and Group Order:",
         "A finite abelian group A has a subgroup of prime index p if and only if its order, |A|, is divisible by p. "
         "Therefore, the abelian group A = G/Φ(G) is filled if and only if its order, |G/Φ(G)|, is divisible by 2 or 3."),
        
        ("6. Relating Orders:",
         "A standard result in group theory is that the set of prime divisors of |G| is the same as the set of prime divisors of |G/Φ(G)|. "
         "This means |G| is divisible by 2 or 3 if and only if |G/Φ(G)| is divisible by 2 or 3."),

        ("7. Conclusion:",
         "By combining these steps, we find that a finite nilpotent group G is filled if and only if its order, |G|, is divisible by the prime 2 or the prime 3.")
    ]

    final_answer = "The finite filled nilpotent groups are precisely the finite nilpotent groups whose order is divisible by 2 or 3."

    # Print the explanation and final answer
    print(f"{title}")
    print("-" * len(title))
    
    for step, text in reasoning:
        print(f"\n{step}\n{text}")

    print("\n" + "="*30)
    print("FINAL ANSWER:")
    print("="*30)
    print(final_answer)


# Execute the function to display the result.
solve_group_theory_problem()