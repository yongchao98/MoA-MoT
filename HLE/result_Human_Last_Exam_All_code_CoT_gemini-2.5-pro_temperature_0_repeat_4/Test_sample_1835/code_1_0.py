def model_generality_constraint():
    """
    This script models the logic of the Generality Constraint, moving from
    understanding a specific proposition 'Fa' to the universal '∀x Fx'.
    """
    # Define a domain of discourse (a set of objects)
    domain_of_discourse = ['a red ball', 'a blue car', 'a red apple', 'a red book']

    # Define a specific object 'a' from the domain
    a = 'a red ball'

    # Define a predicate F(x) which stands for "x is red"
    def F(x):
        return 'red' in x

    # --- Step 1: Understanding 'Fa' ---
    print("--- Step 1: Understanding the proposition 'Fa' ---")
    print(f"Let the predicate F(x) mean 'x is red'.")
    print(f"Let the object 'a' be '{a}'.")
    print("To understand 'Fa' is to evaluate the proposition:")
    print(f"F('{a}') is {F(a)}")
    print("\nThis demonstrates possession of the concept 'F'.")
    print("-" * 50)

    # --- Step 2: Applying the Generality Constraint ---
    print("--- Step 2: Forming the universal proposition '∀x Fx' ---")
    print("Assuming we also understand universal quantification ('For all x...'),")
    print("we can combine our concept 'F' with it to form a new thought.")
    print("We will now evaluate 'For all x in our domain, F(x)'.")
    print(f"Domain = {domain_of_discourse}")

    # --- Step 3: Evaluating '∀x Fx' and printing the final equation ---
    print("\nThe final proposition is an AND-chain of the predicate applied to each object:")

    equation_parts = []
    all_are_true = True

    for item in domain_of_discourse:
        result = F(item)
        # Add each part of the evaluation to our list
        equation_parts.append(f"F('{item}')[{result}]")
        if not result:
            all_are_true = False

    # Join all the parts with 'AND' to form the final equation string
    final_equation = " AND ".join(equation_parts)

    print(f"{final_equation} = {all_are_true}")


if __name__ == '__main__':
    model_generality_constraint()