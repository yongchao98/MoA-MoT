def solve_equivalence_classes():
    """
    Calculates the number of equivalence classes for the given topological space X.

    The space X is the disjoint union of five components. The equivalence relation
    x ~ y holds if there is an auto-homeomorphism of X sending x to y.

    An auto-homeomorphism of X must map connected components to homeomorphic
    connected components. Since all five component spaces of X are topologically
    distinct, any auto-homeomorphism of X maps each component to itself.

    Therefore, the total number of equivalence classes in X is the sum of the
    number of equivalence classes in each component.

    Each of the five component spaces is homogeneous, meaning all its points belong
    to a single equivalence class.
    """

    spaces = [
        {"name": "The torus", "classes": 1},
        {"name": "The sphere", "classes": 1},
        {"name": "The real line", "classes": 1},
        {"name": "A three point discrete space", "classes": 1},
        {"name": "A five point discrete space", "classes": 1},
    ]

    total_classes = 0
    equation_parts = []

    print("The total number of equivalence classes is the sum of the classes from each component.")
    print("Each component is homogeneous, meaning it consists of a single equivalence class.")
    print("-" * 30)

    for space in spaces:
        print(f"Number of classes in '{space['name']}': {space['classes']}")
        total_classes += space['classes']
        equation_parts.append(str(space['classes']))
    
    print("-" * 30)
    print("The final calculation is:")
    
    equation = " + ".join(equation_parts)
    print(f"{equation} = {total_classes}")


if __name__ == "__main__":
    solve_equivalence_classes()